
suppressPackageStartupMessages({
    library(fs)
    library(purrr)
    library(dplyr)
    library(sf)

    library(scater)
    library(splatter)
})

# ##############################################################################
# Define functions for creating grid
# ##############################################################################
make_grid <- function(grid_size, subset=c('full', 'int','ext'))
{
    sfc = sf::st_sfc(sf::st_point(c(0,0)), sf::st_point(c(grid_size,grid_size)))
    grd =  sf::st_make_grid(sfc, n = c(grid_size,grid_size))


    if (subset != "full") {

        border <- get_border(grid_size)
        mask <- if(subset == "int") !border else border

        return(grd[mask])
    }

    grd

}

get_border <- function(grid_size)
{
    border <- which((1:grid_size^2 %% grid_size == 1) | (1:grid_size^2 %% grid_size == 0))
    border <- c(border, 1:grid_size, (grid_size^2-grid_size+1):grid_size^2)


    1:grid_size^2 %in% border

}

## Constants

grid_size <- 120
checker_size <- 15
plot_size <- floor(grid_size/checker_size)
prop <- 0.80

num_genes <- 5000
num_spots <- grid_size^2

# ##############################################################################
# Set up Splatter for simulation of cell and gene matrix
# ##############################################################################
set.seed(13509)

params <- newSplatParams()

# Update number of cells
params <- setParams(
    params,
    batchCells=num_spots,
    nGenes=num_genes,
    group.prob=c(0.5, 0.5) #
)

sim <- splatSimulate(
    params,
    method = "groups",
    verbose = TRUE
)

# Rows - Genes, Columns - Cells
sim_colmeta <- data.frame(colData(sim))
sim_rowmeta <- data.frame(rowData(sim))

sim_colmeta <- sim_colmeta |>
    mutate(Group = ifelse(Group %in% "Group1", "Type1", "Type2")) |>
    group_by(Group) |>
    mutate(
        Group_Order = row_number(),
        CellID = paste0(Group,"_", Group_Order)
    ) |>
    ungroup() |>
    mutate(
        Matrix_Order = row_number()
    )

# Used to generate grid
group_tally <- sim_colmeta %>% group_by(Group) %>% tally()
n_type1 <- group_tally[[1,2]]

dimnames(sim)[[2]] <- sim_colmeta$CellID

# ##############################################################################
# Create grid and merge metadata
# ##############################################################################

sfc <- st_sf(geometry=make_grid(grid_size, 'full'))

sfc$spatial_order <- seq(grid_size^2)
sfc$row <- floor((sfc$spatial_order - 1)/grid_size)
sfc$col <- (sfc$spatial_order-1) %% grid_size

checker_inds <- split(seq(0,grid_size-1),sort(rep(seq(plot_size),checker_size)))
checker <- unlist(checker_inds[seq(plot_size-1,by=-2)])

sfc <- sfc %>%
    mutate(Checker = case_when(
        row %in% checker & col %in% checker ~ 'Black',
        (!row %in% checker) & (!col %in% checker) ~ 'Black',
        .default = 'White'
    ))

# --proportion of Checker 1 will be Type 1
# Use `n=` argument to match dimensions of the count matrix
set.seed(13509)
n_type1_prop = floor(prop*n_type1)
n_type1_min = n_type1-n_type1_prop

inds_type1 <- sfc |>
    ungroup() |>
    group_split(Checker) |>
    map2_dfr(c(n_type1_prop, n_type1_min), ~ slice_sample(.x, n = .y)) %>%
    pull(spatial_order)

sfc <- sfc %>%
    mutate(Group = ifelse(spatial_order %in% inds_type1,"Type1","Type2"))

# Bind coordinate to CellID
sfc <- sfc %>%
    group_by(Group) |>
    mutate(
        Group_Order_sp = row_number(),
        CellID = paste0(Group,"_", Group_Order_sp)
    ) |>
    data.frame() |>
    rename_with(~ paste0(.x, "_", prop), ends_with("spatial_order")) |>
    rename_with(~ paste0(.x, "_", grid_size, "_", checker_size), ends_with("Checker")) |>
    select(-geometry, CellID, contains("spatial_order"),
           contains("Checker"), Group, row, col)

# coords <- sfc |> select(contains("spatial_order"), contains("Checker"), row, col)

cell_metadata <-
    merge(sim_colmeta, sfc, by=c('CellID','Group')) |>
    arrange(Matrix_Order) %>%
    select(-Batch, -contains("Group_Order"))

# Save metadata
outfile <- path_wd('checker/cell_metadata.csv')
write.csv(cell_metadata, file=outfile, row.names=F)
