using Nucleus

using Omics

using Human

# ---- #

const DI = joinpath(dirname(Human.PK), "Data.pr", "output", "gene_set")

const D1 = Nucleus.Dictionary.rea(joinpath(DI, "immune_population", "all.json"))

const D2 = Nucleus.Dictionary.rea(joinpath(DI, "neural_population", "all.json"))

#pop!(D1, "Progenitors")

# ---- #

const ST_ = intersect(
    readlines(joinpath(Human.OU, "1.tsv")),
    vcat(reduce(vcat, values(D1)), reduce(vcat, values(D2))),
)

# ---- #

const D3 = Nucleus.Table.make_dictionary(
    Nucleus.Table.rea(Omics.Gene.HG),
    ["symbol"],
    "gene_group",
)

# ---- #

const D4 = Dict(
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

# ---- #

function text(di, s1)

    s2_ = String[]

    for (s2, s1_) in di

        if s1 in s1_

            push!(s2_, s2)

        end

    end

    Nucleus.Strin.text(s2_)

end

Nucleus.Table.writ(
    joinpath(Human.OU, "2.tsv"),
    Nucleus.Table.make(
        "Gene",
        ST_,
        ["Group", "Immune Population", "Neural Population", "Top Blastp against Mouse"],
        hcat(
            map(st -> get(D3, st, ""), ST_),
            map(st -> text(D1, st), ST_),
            map(st -> text(D2, st), ST_),
            map(st -> get(D4, st, NaN), ST_),
        ),
    ),
)
