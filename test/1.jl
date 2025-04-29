using Nucleus

using Omics

using Human

# ---- #

const A = Nucleus.Table.rea(Omics.Gene.EN)

write(
    joinpath(Human.OU, "1.tsv"),
    join(
        sort!(
            setdiff(
                skipmissing(
                    A[findall(==("protein_coding"), A[!, "Gene type"]), "Gene name"],
                ),
                skipmissing(
                    Nucleus.Table.rea(joinpath(Human.IN, "ensembl.mouse_human.tsv.gz"))[
                        !,
                        "Human gene name",
                    ],
                ),
            ),
        ),
        '\n',
    ),
)
