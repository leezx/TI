# see: https://github.com/dynverse/ti_dpt
# see: https://github.com/dynverse/dynwrap

# get my data
# library(Seurat)
# print(load("annotated.E95_P21.seurat.Rdata"))
# for (i in c("E9.5","E10.5","E11.5","E12.5","E13.5","E16.5","P19","P21")) {
#     # print(i)
#     tmp.stage <- i
#     tmp.stage <- "E13.5"
#     print(tmp.stage)
#     tmp.cells <- names(ens.integrated$seurat_clusters)[ens.integrated$stage == tmp.stage]
#     #
#     set.seed(49)
#     if (length(tmp.cells) > 3000) {
#         tmp.cells <- sample(tmp.cells, size = 3000)
#     }
#     print(length(tmp.cells))
#     #
#     tmp.counts <- ens.integrated@assays$RNA@counts[,tmp.cells]
#     tmp.expression <- ens.integrated@assays$RNA@scale.data[,tmp.cells]
#     tmp.subtype <- ens.integrated[,tmp.cells]$lineage.sub
#     #
#     # save(tmp.counts, tmp.expression, tmp.subtype, file=paste("stages/",i,"_",length(tmp.cells),"cells.Rdata",sep=""))
#     break
# }

my_ti_dpt <- function(tmp.expression) {
    
    print(dim(tmp.expression))

    # connect to API
    library(dplyr, warn.conflicts = FALSE)
    library(purrr, warn.conflicts = FALSE)
    library(destiny, warn.conflicts = FALSE)

    expression <- t(tmp.expression) %>% as.matrix

    # run diffusion maps
    dm <- destiny::DiffusionMap(
      data = expression,
      sigma = "local",
      distance = "euclidean",
      n_eigs = 20,
      density_norm = T,
      n_local = c(5,7)
      #vars = parameters$features_id
    )

    tips <- random_root(dm)

    dpt <- destiny::DPT(
      dm,
      w_width = 0.1,
      tips = tips
    )

    # find DPT tips
    tips <- destiny::tips(dpt)
    tip_names <- rownames(expression)[tips]

    #####################################
    ###     SAVE OUTPUT TRAJECTORY    ###
    #####################################
    cell_ids <- rownames(expression)

    # construct grouping
    grouping <- dpt@branch[,1] %>%
      ifelse(is.na(.), 0, .) %>%
      as.character() %>%
      paste0("Tip", .) %>%
      purrr::set_names(cell_ids)

    group_ids <- sort(unique(grouping))

    # collect distances from tips
    tip_dists <- dpt[,tips] %>%
      magrittr::set_colnames(., paste0("Tip", seq_len(ncol(.)))) %>%
      magrittr::set_rownames(cell_ids)

    # calculate progressions
    outs <- map(
      group_ids,
      function(gid) {
        cat("Processing ", gid, "\n", sep = "")
        cixs <- which(grouping == gid)
        cids <- cell_ids[cixs]
        if (length(cids) > 0) {
          if (gid == "Tip0") {
            progr <- data_frame(
              cell_id = cids,
              from = gid,
              to = sample(setdiff(group_ids, gid), length(cids), replace = TRUE),
              percentage = 0
            )
            list(progr = progr, milnet = NULL)
          } else {
            # calculate min dist of gid to all other cells
            max_range <- min(tip_dists[-cixs, gid])

            # calculate percentage value of cells in gid
            percentage <- 1 - pmin(tip_dists[cixs, gid] / max_range, 1)
            progr <- data_frame(cell_id = cids, from = "Tip0", to = gid, percentage = percentage)
            milnet <- data_frame(from = "Tip0", to = gid, length = max_range, directed = FALSE)
            list(progr = progr, milnet = milnet)
          }
        } else {
          list()
        }
      }
    )
    progressions <- map_df(outs, ~ .$progr)
    milestone_network <- map_df(outs, ~ .$milnet)

    # collect dimred
    dimred <- dm@eigenvectors %>%
      magrittr::set_colnames(., paste0("Comp", seq_len(ncol(.)))) %>%
      magrittr::set_rownames(cell_ids)

    library(dynwrap)
    output <-
      wrap_data(
        cell_ids = cell_ids
      ) %>%
      add_grouping(
        group_ids = group_ids,
        grouping = grouping
      ) %>%
      add_trajectory(
        milestone_ids = group_ids,
        milestone_network = milestone_network,
        progressions = progressions
      ) %>% 
      add_dimred(
        dimred = dimred
      ) %>%
      add_timings(
        timings = timings
      )


    # # connect to dyno
    # model <- output 
    # model <- model %>% add_dimred(dyndimred::dimred_mds, expression_source = expression)
    # library(ggplot2)

    # options(repr.plot.width=8, repr.plot.height=6)
    # plot_dimred(
    #   model,  
    #   expression_source = expression, 
    #   grouping = tmp.subtype
    # ) + ggtitle("DPT")
    # ## Coloring by grouping

    return(output)
}    
    
# ggsave(filename = paste("stages/E13.5_DPT.pdf"), width = 8, height = 6, dpi = 800)
