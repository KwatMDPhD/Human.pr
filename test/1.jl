using Nucleus

using Omics

using Human

# ---- #

const BO = false

const A1 = Nucleus.Table.rea(Omics.Gene.EN)

const A2 = Nucleus.Table.rea(joinpath(Human.IN, "ensembl.mouse_human.tsv.gz"))

const ST_ = sort!(
    setdiff(
        skipmissing(A1[findall(==("protein_coding"), A1[!, "Gene type"]), "Gene name"]),
        skipmissing(A2[!, "Human gene name"]),
    ),
)

write(joinpath(Human.OU, "human_gene.tsv"), join(ST_, '\n'))
