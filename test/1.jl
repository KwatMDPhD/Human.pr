using Nucleus

using Omics

using Human

# ---- #

const S1_, S2_ = eachcol(Nucleus.Table.rea(Omics.Gene.EN)[!, ["Gene name", "Gene type"]])

write(
    joinpath(Human.OU, "1.tsv"),
    join(
        sort!(
            setdiff(
                skipmissing(S1_[findall(==("protein_coding"), S2_)]),
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
