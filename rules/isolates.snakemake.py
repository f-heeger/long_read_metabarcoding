rule collectIsoData:
    input: "all_otu{ident}_table.tsv"
    output: tab="isolates/isolate{ident}_table.tsv"
    script:
        "../scripts/isolates_collectIsoData.py"

rule plotIsoData:
    input: "isolates/isolate{ident}_table.tsv"
    output: ssu="isolates/isolate{ident}_ssu.svg", its="isolates/isolate{ident}_its.svg", lsu="isolates/isolate{ident}_lsu.svg"
    conda:
        "../envs/ggplot.yaml"
    script:
        "../scripts/isolates_plotIsoClassification.R"
