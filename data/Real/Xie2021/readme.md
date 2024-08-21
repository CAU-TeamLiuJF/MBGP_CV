
# Data Source

[Xie, L.](https://doi.org/10.1111/age.13121), J. Qin, L. Rao, X. Tang and D. Cui et al., 2021 Accurate prediction and genome-wide association analysis of digital intramuscular fat content in longissimus muscle of pigs. Animal Genetics 52: 633-644.

[source URL](https://doi.org/10.6084/m9.figshare.14274785.v1)

# Data Overview

## Populations

A total of **1709** pigs were randomly selected from **multiple pig populations**, including **228** Landrace pigs (141 females and 87 castrates), **641** Large White pigs (YK, **407** females and **234** castrates), **605** Landrace × Large White crossbred pigs (LY, 355 females and 250 castrates), and **235** Duroc × Landrace × Large White hybrid pigs (DLY, 128 females and 107 castrates). All pigs were reared at **Muyuan** company, which operates under a large-scale intensive farming system. These pigs were raised in a consistent rearing environment with free access to the same commercial feed and water. Diets were formulated according to the recommendations of the **Chinese National Feed Standard** (GB 2004), tailored to meet the needs of pigs at different growth stages.

| Breed | Number of Individuals |
| ---- | ------ |
| YY   | 641    |
| LL   | 228    |
| LY   | 601    |
| DLY  | 232    |

## Phenotypes

All pigs were fed a growing diet from 41 to 120 days of age, followed by a finishing diet from 121 to 180 days of age (Table S1). At 180 days, the pigs were **slaughtered uniformly** at Longda Muyuan Meat Co., Ltd. (Henan, China). In short, following the protocol described in the "Livestock and Poultry Slaughter Operation Procedures - Pigs" (GB/T 17236-2019), pigs were processed through CO2 stunning, bleeding, cleaning, head and foot removal, dehairing, evisceration, cutting, weighing, blast freezing, and storage. After blast cooling at temperatures between −10°C and −15°C for 2–4 hours, the surface temperature of the carcass was close to −2°C. The carcass was then transferred to a 4°C cold storage room and cooled for at least 20 hours. After 24 hours of chilling, the longissimus dorsi (LD) muscle at the fifth to sixth rib (5-7 cm thick) of the left half of the carcass was sampled. Finally, 2.5 cm LD slices, with the subcutaneous and connective tissue removed, were used for computer **image** acquisition, routine information processing, and **feature extraction**.

After trimming the subcutaneous and connective tissue from the pork slice, digital image processing technology was used to collect images via a **homemade image acquisition system**. The constructed image acquisition system consisted of four parts (Figure 1): a black box (40 cm × 40 cm × 40 cm) with reflective polyethylene material placed around it to help distribute light evenly; a circular fluorescent lamp (Philips T5-40W/6500K; Royal Dutch Philips Electronics Ltd.); a color digital camera (model D5200; Nikon Corporation); and a laptop (ACER E5-572G, CPU: Intel Core i5 4210M, Graphics: NVIDIA GeForce 840M, RAM: 8G). The LD slicer was then placed in the center of the black box (40 cm below the camera lens) for image acquisition (Figure 2a: sample image). Images were processed and analyzed using Matlab software (version R2016a; The MathWorks).

| MS   | Number of samples | PFAI (%)      |        |      |       | P value    |
| ---- | ----------------- | ------------- | ------ | ---- | ----- | ---------- |
|      |                   | Mean ± SD     | Median | Min  | Max   |            |
| 1    | 565               | 0.10 ± 0.41a  | 1.02   | 0    | 1.99  | <2 × 10−16 |
| 1.5  | 688               | 2.13 ± 0.48b  | 2.08   | 1.31 | 3.46  |            |
| 2    | 321               | 3.41 ± 0.41c  | 3.43   | 2.7  | 4.46  |            |
| 2.5  | 59                | 4.73 ± 0.41d  | 4.66   | 4.01 | 5.49  |            |
| 3    | 47                | 5.91 ± 0.41e  | 5.8    | 5.28 | 6.88  |            |
| 3.5  | 18                | 7.40 ± 0.56f  | 7.22   | 6.56 | 8.37  |            |
| 4    | 9                 | 9.89 ± 1.15g  | 9.82   | 8.07 | 11.75 |            |
| 4.5  | 2                 | 13.48 ± 1.82h | 13.48  | 12.2 | 14.77 |            |

- PFAI, proportion of fat area in the image; MS, marbling score; SD, standard deviation.
- Within each row, different letters (a–h) denote statistical differences between groups (*P* < 0.05). *P* value: highly significant after multiple-group comparisons.

## Genotypes

Genomic DNA was extracted from the **muscle** tissue of each animal using a standard phenol/chloroform extraction method and dissolved in Tris-EDTA buffer. All DNA samples were quantified and diluted to a final concentration of 75 ng/μl. Genotyping was performed using the CC1 **PorcineSNP50 BeadChip** (containing **51,368** SNPs) following the manufacturer's protocol (Sscrofa10.2).

Quality control was performed to exclude SNPs with a **call rate** <95%, **MAF** <5%, or Hardy-Weinberg equilibrium P-value <10^-5. The call rate for all animals was >90%. A total of **40,016** SNPs and **1709** animals passed quality control and were retained for further statistical analysis.
