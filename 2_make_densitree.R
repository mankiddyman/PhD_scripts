library(phangorn)
library(ape)


tree_dir <- "/netscratch/dep_mercier/grp_marques/Aaryan/Drosera/genespace/reg_orthogroups_phasing/orthofinder/Results_Dec23/Single_Copy_Orthogroups_Trees/Single_Copy_Orthogroups_Gene_Trees/Sanitized_Names"

tree_files <- list.files(
  tree_dir,
  pattern = "\\.txt$|\\.nwk$|\\.tree$",
  full.names = TRUE
)

trees <-  lapply(
  tree_files,
  read.tree
)
# Define Droseraceae taxa (must match tip labels exactly)
droseraceae <- c(
  "D_regia",
  "Dio_muscipula",
  "D_capensis"
  )

class(trees) <- "multiPhylo"
class(trees)      # "multiPhylo"
length(trees)     # ~400

# Root all trees so that the basal split separates Droseraceae
trees <- lapply(seq_along(trees), function(i) {

  tr <- trees[[i]]

  present <- droseraceae[droseraceae %in% tr$tip.label]

  if (length(present) < 2) {
    stop(
      paste0("Tree ", tree_files[i],
             " has fewer than 2 Droseraceae taxa — cannot root")
    )
  }

  mrca_node <- getMRCA(tr, present)

  if (is.null(mrca_node)) {
    stop(
      paste0("Droseraceae not monophyletic in tree: ",
             tree_files[i])
    )
  }

  root(tr, node = mrca_node, resolve.root = TRUE)
})

class(trees) <- "multiPhylo"

# santity check: all trees should have the same set of tips 
ref_tips <- sort(trees[[1]]$tip.label)


stopifnot(

  all(sapply(trees, function(t)

    identical(sort(t$tip.label), ref_tips)))

)



# sanity check: all trees should have branch lengths
stopifnot(
  all(sapply(trees, function(t) !is.null(t$edge.length)))
)

# MAKE A FOLDER CALLED ROOTED_TREES AND SAVE THE ROOTED TREES THERE
rooted_tree_dir <- file.path(tree_dir, "ROOTED_TREES")
if (!dir.exists(rooted_tree_dir)) {
  dir.create(rooted_tree_dir)
}
write.tree(trees, file = file.path(rooted_tree_dir, "rooted_trees.nwk"))

# Plot the DensiTree
pdf("Drosera_orthogroups_densitree.pdf", width = 10, height = 10)
densiTree(trees,type='cladogram',width=2,alpha=0.3)
dev.off()
