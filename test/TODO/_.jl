pa1 = "set1_genes.json"

#pa2 = "set2_genes.json"

na1 = "Immune Population"

#na2 = "Neural Population"

na1c = Omics.path.clean(na1)

#na2c = Omics.path.clean(na2)

se1_ge_ = Omics.dict.read(joinpath("..", "input", pa1))

pop!(se1_ge_, "Progenitors")

#se2_ge_ = Omics.dict.read(joinpath("..", "input", pa2))

ge1_ = vcat(values(se1_ge_)...)

#ge2_ = vcat(values(se2_ge_)...)

length(ge1_)

#length(ge2_)

# ---- #

na = na1

nac = na1c

se_ge_ = se1_ge_

ge_ = ge1_

length(ge_)

# ---- #

if BO

    na = "na2"

    nac = Omics.path.clean(na)

    se_ge_ = sort(merge(se1_ge_, se2_ge_))

    ge_ = vcat(values(se_ge_)...)

    println(length(ge_))

end

# ---- #

di = joinpath("..", "output", nac)

mkpath(di)

# ---- #

co_ = intersect(hu_, ge_)

# ---- #

se_ = [join([se for (se, ge_) in se_ge_ if co in ge_], " & ") for co in co_]

# ---- #

hg = Omics.gene.read_hgnc()

ge_gr = Omics.data_frame.map_to_column(hg, ["symbol", "gene_group"])

ge_na = Omics.data_frame.map_to_column(hg, ["symbol", "name"])

gr_ = [get(ge_gr, co, uppercasefirst(ge_na[co])) for co in co_]

# ---- #

ge_bl = Dict(
    "GNLY" => 0.292,
    "COL9A3" => 0.512,
    "FAM177B" => 0.388,
    "GLYATL1B" => 0.412,
    "ZNF860" => 0.438,
    "PLEKHG7" => 0.258,
    "GCSAML" => 0.292,
    "IL3RA" => 0.305,
    "KRT81" => 0.881,
    "KRT86" => 0.872,
    "OR14L1P" => 0.530,
    "OR6K3" => 0.784,
    "OXER1" => 0.418,
    "RAB32" => 0.839,
    "TRIM49D1" => 0.364,
    "TRIM49D2" => 0.364,
    "TRIM51" => 0.367,
    "TRIM64" => 0.348,
    "TRIM64B" => 0.348,
    "UGT2B11" => 0.704,
    "ADGRE3" => 0.459,
    "CFH" => 0.614,
    "CXCL2" => 0.735,
    "CXCL3" => 0.685,
    "GPR84" => 0.856,
    "VSTM1" => 0.382,
    "NT5DC4" => 0.654,
    "CASP5" => 0.636,
    "CTSL" => 0.785,
    "LILRA1" => 0.507,
    "LILRB2" => 0.452,
    "LILRA2" => 0.513,
    "PILRA" => 0.469,
    "SIGLEC14" => 0.498,
    "TICAM2" => 0.744,
    "KIR2DP1" => 0.401,
    "KLRF1" => 0.303,
    "LGALS9B" => 0.685,
    "LGALS9C" => 0.688,
    "APOBEC3A" => 0.426,
    "BTNL8" => 0.396,
    "CEACAM3" => 0.528,
    "CEACAM4" => 0.472,
    "CSNK1A1L" => 0.910,
    "CXCL1" => 0.701,
    "CXCR1" => 0.686,
    "FCAR" => 0.327,
    "H2AC19" => 1,
    "H3C4" => 1,
    "HSPA6" => 0.815,
    "LINC02218" => 0.441,
    "PEAK3" => 0.268,
    "PI3" => 0.941,
    "S100A12" => 0.388,
    "S100P" => 0.545,
    "SIGLEC5" => 0.455,
    "SIRPB2" => 0.348,
    "SIRPD" => 0.48,
    "ST20" => 0.724,
    "TMEM272" => 0.368,
    "APOBEC3B" => 0.359,
    "PLAAT2" => 0.701,
    "ZNF683" => 0.560,
    "TRABD2A" => 0.648,
    "FAAH2" => 0.304,
    "CD1B" => 0.489,
    "CD1E" => 0.521,
    "CLIC2" => 0.662,
    "FPR3" => 0.634,
    "CCDC183" => 0.76,
    "CTSV" => 0.792,
    "GGTA1" => 0.654,
    "KCNK17" => 0.42,
    "NLRP7" => 0.408,
    "NOTCH4" => 0.808,
    "PROC" => 0.695,
    "SLC9C2" => 0.292,
    "SMIM6" => 0.565,
)

bl_ = [get(ge_bl, co, NaN) for co in co_]

# ---- #

ge_x_io = sort(
    DataFrame(
        "Gene" => co_,
        uppercasefirst(na) => se_,
        "Group" => gr_,
        "Top blastp against mouse" => bl_,
    ),
    na,
)

if BO

    ge_x_io = ge_x_io[
        [
            !any(all(haskey(sen_ge_, sp) for sp in sp_) for sen_ge_ in (se1_ge_, se2_ge_)) for sp_ in split.(ge_x_io[!, na], " & ")
        ],
        :,
    ]

end

Omics.table.write(joinpath(di, "gene_x_information.tsv"), ge_x_io)

ge_x_io

# ---- #

he = 1500

if BO

    he /= 2

end

wi = round(Int, he * MathConstants.golden)

axis_title_font_size = 48

layout = Dict(
    "barmode" => "stack",
    "height" => he,
    "width" => wi,
    "margin_b" => 240,
    "title" => Dict("text" => "Number of Human Genes in $(na)", "font_size" => 32),
    "yaxis" => Dict(
        "title" => Dict(
            "text" => "Number of Human Genes",
            "font_size" => axis_title_font_size,
        ),
        "automargin" => true,
    ),
    "xaxis" => Dict(
        "title" => Dict("text" => "$na", "font_size" => axis_title_font_size),
        #"tickangle" => 60,
        "automargin" => true,
    ),
)

# ---- #

se_n = sort(countmap(ge_x_io[!, na]); byvalue = true)

x = collect(keys(se_n))

y = collect(values(se_n))

Omics.figure.plot_bar([y], [x]; layout = layout, ou = joinpath(di, "genes.html"))

# ---- #

name_ = []

y_ = []

id = 1

for gr in unique(ge_x_io[!, "Group"])

    push!(name_, gr)

    n_be_ = []

    for se in x

        da = ge_x_io[ge_x_io[!, na] .== se, :]

        gr_ = da[!, "Group"]

        n_be = sum(gr_ .== gr)

        push!(n_be_, n_be)

    end

    push!(y_, n_be_)

end

Omics.figure.plot_bar(
    y_,
    fill(x, length(y_));
    name_ = name_,
    layout = layout,
    ou = joinpath(di, "genes_by_group.html"),
)

# ---- #

if BO

    se2_se1_ = Dict()

    for st in ge_x_io[!, na]

        sp_ = split(st, " & ")

        se2_ = [haskey(se2_ge_, sp) for sp in sp_]

        no_ = sp_[.!se2_]

        for se2 in sp_[se2_]

            if haskey(se2_se1_, se2)

                append!(se2_se1_[se2], no_)

            else

                se2_se1_[se2] = copy(no_)

            end

        end

    end

    se2_se1_ = sort(se2_se1_)

    la["xaxis"]["title"]["text"] = "$na1"

    x = collect(keys(sort(countmap(vcat(values(se2_se1_)...)); byvalue = true)))

    Omics.figure.plot_bar(
        [[sum(se1_ .== xx) for xx in x] for (se2, se1_) in se2_se1_],
        fill(x, length(se2_se1_));
        name_ = collect(keys(se2_se1_)),
        layout = layout,
        ou = joinpath(di, "genes_by_set.html"),
    )

end
