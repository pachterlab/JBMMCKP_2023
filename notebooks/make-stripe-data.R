# Install Splatter
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
# BiocManager::install("splatter")

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
checker_size <- 30
plot_size <- floor(grid_size/checker_size)


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

sfc <- sfc |> arrange(row)

sfc["stripe"] <- sort(rep(seq(plot_size), checker_size*grid_size))


# ##############################################################################
# Assign cell types to stripes
# ##############################################################################
seqs <- c(0.25, 0.40, 0.60,0.75)

# Assign this many locs `Type1` in each stripe
# If the fractions aren't pretty, add them to the last stripe
seqs_n <- floor(checker_size*grid_size*seqs)

differ <- n_type1-sum(seqs_n)
seqs_n[plot_size] <- seqs_n[plot_size] + differ

# Get the spatial indicies of `Type1` locations
inds_type1 <- sfc |>
    ungroup() |>
    group_split(stripe) |>
    map2_dfr(seqs_n, ~ slice_sample(.x, n = .y)) |>
    pull(spatial_order)

sfc <- sfc %>%
    mutate(Group = ifelse(spatial_order %in% inds_type1,"Type1","Type2"))
colors <- c("Type2" = "#E0AFCA","Type1"="#B14380")

# plot(sfc$row, sfc$col, col=colors[sfc$Group])

# Bind Coordinate to CellID
sfc <- sfc |>
    group_by(Group) |>
    mutate(
        Group_Order_sp = row_number(),
        CellID = paste0(Group,"_", Group_Order_sp)
    ) |>
    data.frame() |>
    rename_with(~ paste0(.x, "_", grid_size, "_", checker_size), ends_with("stripe")) |>
    select(-geometry, CellID, contains("spatial_order"),
           contains("stripe"), Group, row, col)

cell_metadata <-
    merge(sim_colmeta, sfc, by=c('CellID','Group')) |>
    arrange(Matrix_Order) %>%
    select(-Batch, -contains("Group_Order"))


# Save metadata
outfile <- path_wd('stripe/cell_metadata.csv')
write.csv(cell_metadata, file=outfile, row.names=F)
