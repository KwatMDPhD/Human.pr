using Nucleus

using Human

# ---- #

const _, _, ST_, S =
    Nucleus.Table.ge(Nucleus.Table.rea(joinpath(Human.OU, "2.tsv"))[!, 1:3])

# ---- #

function make(st)

    "title" => Dict("text" => st)

end

for nd in eachindex(ST_)

    st = ST_[nd]

    st_ = skipmissing(S[:, nd])

    un_ = unique(st_)

    um_ = map(un -> count(==(un), st_), un_)

    in_ = sortperm(um_)

    Nucleus.Plotly.writ(
        joinpath(Human.OU, "$st.html"),
        (Dict("type" => "bar", "y" => um_[in_], "x" => un_[in_]),),
        Dict(
            "barmode" => "stack",
            make("Human Genes"),
            "yaxis" => Dict(make("Count")),
            "xaxis" => Dict(make(st)),
        ),
    )

end
