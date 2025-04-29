using Nucleus

using Human

# ---- #

const A = Nucleus.Table.rea(joinpath(Human.OU, "2.tsv"))

# ---- #

for st in ("Group", "Immune Population", "Neural Population")

    st_ = skipmissing(A[!, st])

    un_ = unique(st_)

    um_ = map(un -> count(==(un), st_), un_)

    in_ = sortperm(um_)

    Nucleus.Plotly.writ(
        joinpath(Human.OU, "$st.html"),
        (Dict("type" => "bar", "y" => um_[in_], "x" => un_[in_]),),
        Dict(
            "width" => max(1600, lastindex(un_) * 16),
            "barmode" => "stack",
            "title" => Dict("text" => "Human Genes"),
            "yaxis" => Dict("title" => Dict("text" => "Count")),
            "xaxis" => Dict("title" => Dict("text" => st)),
        ),
    )

end
