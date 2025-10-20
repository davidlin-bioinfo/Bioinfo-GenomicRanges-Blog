# Practical TSS Variant Analysis Using GenomicRanges: A Focus on Dual Promoters in Cancer

May 4, 2025

Zhengyi Lin

<aside>
🔖 **In this article, we’ll walk through a practical case study on transcription start site (TSS) analysis to explore how to use the R package GenomicRanges for assessing pathogenicity in dual promoters and their upstream/downstream regions.**

</aside>

## `*Background`(72 Words)*

In cancer research, variations at transcription start sites (TSSs) play a crucial role in regulating gene expression. Recent studies show the presence of dual promoters in cancer cells—meaning a single gene may initiate transcription from two distinct TSSs[Fig1][1]. Compared to the classic single-promoter model, dual promoters can give rise to alternative mRNA transcripts, potentially producing functionally diverse proteins that may influence cancer cell proliferation and behaviour[2].

Our goal in this analysis is to use the GenomicRanges package to investigate variations in these dual TSS regions, with a particular focus on their overlap with known pathogenic variants. Leveraging the tools within GenomicRanges, we analyse how these dual promoter-associated TSSs intersect with variants reported in public datasets such as ClinVar and 100k Genomes Project. Ultimately, we aim to assess the potential roles these variants may play in cancer.

![Fig1. Different types of transcription initiation patterns in gene promoters[1]](Practical%20TSS%20Variant%20Analysis%20Using%20GenomicRanges%201e836f2bf90a80b89430fc6548ed29f2/Screenshot_2025-05-03_at_07.59.51.png)

Fig1. Different types of transcription initiation patterns in gene promoters[1]

## `*So, What Exactly Is GenomicRanges?（186 Words）*`

Let’s take a closer look at the R package *GenomicRanges*, a powerful tool designed for working with genomic interval data. Widely used in bioinformatics, it’s ideal for tasks like variant analysis, transcriptomic data processing, and ChIP-seq region handling[3].
At its core, GenomicRanges helps efficiently manage genomic intervals—such as chromosome coordinates, start/end positions, strand orientation, and metadata.

**First, Two Key Concepts**

1️⃣ Genomic Region
A genomic region refers to a continuous segment on a chromosome, typically defined by:
• Chromosome name
• Start and end positions
• Strand (+/–)

In our project, we filtered data based on experimental needs—excluding sex chromosomes, keeping regions under 150 bp, and focusing on segments with dual transcription start sites. This gave us a clean CSV file for further analysis.

![Fig. 2: The input CSV file after filtering](Practical%20TSS%20Variant%20Analysis%20Using%20GenomicRanges%201e836f2bf90a80b89430fc6548ed29f2/Screenshot_2025-05-02_at_16.54.23.png)

Fig. 2: The input CSV file after filtering

2️⃣ **GRanges Object**

This is the heart of GenomicRanges. Think of a GRanges object as a “smart” version of a VCF file—structured, flexible, and easy to manipulate.

A key step is to clearly define the information needed for analysis—for example, after loading the CSV, we created a GRanges object using coordinates, strand, gene names, and a YC/YR promoter dominance label to enable downstream stratification.

![Screenshot 2025-05-06 at 15.17.39.png](Practical%20TSS%20Variant%20Analysis%20Using%20GenomicRanges%201e836f2bf90a80b89430fc6548ed29f2/9b39c24b-6122-4ca8-a97d-373aa840e977.png)

## **`*Applying GenomicRanges to Analysis`**(459 Words)*

We will first define the regions, then perform the overlap analysis, and finally visualize the results through plotting.

### **Defining Genomic Regions**

When analyzing TSS regions, we not only focus on the TSS itself but also consider the surrounding upstream and downstream regions, as they may contain key regulatory elements. To comprehensively study transcriptional regulation, we expand the TSS regions to include 100 base pairs upstream and downstream as flanking sites. As a result, we define two comparison groups: dual-promoter TSS regions and flanking sites.

![Screenshot 2025-05-06 at 15.18.51.png](Practical%20TSS%20Variant%20Analysis%20Using%20GenomicRanges%201e836f2bf90a80b89430fc6548ed29f2/Screenshot_2025-05-06_at_15.18.51.png)

In this code, the `flank()` function extends the TSS regions by 100 base pairs both upstream and downstream. By setting `ignore.strand = FALSE`, we ensure strand orientation is respected, which is essential for correctly defining the upstream and downstream regions depending on whether the gene is on the forward or reverse strand.
You can also apply other GenomicRanges functions—such as `reduce()` for merging regions, `subset()` for filtering, `distance()` for measuring distances between intervals to further refine and explore your results.

### **Overlap Analysis**

**GenomicRanges** excels at efficiently computing overlaps between genomic intervals, which is essential in many genomic workflows. In this case, we identify variants from a **VCF file** (e.g., ClinVar) that overlap with or are near our TSS regions to explore their potential impact on gene transcription.

Here I using TSS positions as examples:

![Screenshot 2025-05-06 at 15.19.15.png](Practical%20TSS%20Variant%20Analysis%20Using%20GenomicRanges%201e836f2bf90a80b89430fc6548ed29f2/Screenshot_2025-05-06_at_15.19.15.png)

One of the major advantages of using **GenomicRanges** is that it shares a unified data structure with other **Bioconductor** packages(e.g. VariantAnnotation). This allows for seamless conversion of your VCF file into a **GRanges** object once loaded, enabling easy downstream analysis. The `findOverlaps()` function identifies overlapping regions between TSS and variant positions, and the returned index pairs let us quickly extract the relevant entries from both datasets for further functional or clinical investigation.

### **Annotation**

Once we’ve find overlaps, the next step is to annotate these variants with biologically meaningful information—such as clinical significance. 

![Screenshot 2025-05-06 at 15.19.52.png](Practical%20TSS%20Variant%20Analysis%20Using%20GenomicRanges%201e836f2bf90a80b89430fc6548ed29f2/Screenshot_2025-05-06_at_15.19.52.png)

We using data from **ClinVar**, then extract annotations CLNSIG from info column, which indicates whether a variant has been linked to disease.

### **Visualising the Results**

To make the findings more intuitive, we can visualise the distribution of variants—broken down by promoter type (TSS vs Flanking sites) and annotated significance. 

![Screenshot 2025-05-06 at 15.20.37.png](Practical%20TSS%20Variant%20Analysis%20Using%20GenomicRanges%201e836f2bf90a80b89430fc6548ed29f2/Screenshot_2025-05-06_at_15.20.37.png)

Through this approach, we not only pinpoint variants within TSS regions, but also assess their disease relevance and distribution across promoter categories—offering a more complete picture of their potential functional roles. Proportion of pathogenic variants is markedly higher at dual-promoter TSS regions (27.9%) compared to flanking regions (9.5%), suggesting stronger disease association at promoter sites.[Fig.2]

![Fig. 2](Practical%20TSS%20Variant%20Analysis%20Using%20GenomicRanges%201e836f2bf90a80b89430fc6548ed29f2/Screenshot_2025-05-03_at_11.04.17.png)

Fig. 2

![Screenshot 2025-05-03 at 11.04.10.png](Practical%20TSS%20Variant%20Analysis%20Using%20GenomicRanges%201e836f2bf90a80b89430fc6548ed29f2/Screenshot_2025-05-03_at_11.04.10.png)

### `*In the end`(41 Words)*

At this stage, we have successfully performed chromosomal region comparisons using the *GenomicRanges* R package and obtained preliminary experimental results. Building on this analysis framework, we will next extend our investigation from germline cells to somatic cells for further in-depth analysis.

### References

1. Nepal C, Hadzhiev Y, Balwierz P, Tarifeño-Saldivia E, Cardenas R, Wragg JW, et al. Dual-initiation promoters with intertwined canonical and TCT/TOP transcription start sites diversify transcript processing. Nat Commun. 2020 Jan 10;11(1):168.
2. Wragg JW, White PL, Hadzhiev Y, Wanigasooriya K, Stodolna A, Tee L, et al. Intra-promoter switch of transcription initiation sites in proliferation signaling-dependent RNA metabolism. Nature Structural & Molecular Biology. 2023 Nov 23;30(12):1970.
3. Lawrence M, Aboyoun P, Pages H, and Gentleman R (2013). *GenomicRanges: Representation and manipulation of genomic intervals*. R package version 1.61.0 https://bioconductor.org/packages/GenomicRanges
