Code for upcoming paper, titled:

# The impact of mammal extinctions on pathogen richness in extant hosts.

by Tomos O. Prys-Jones<sup>1</sup>, Andrew J. Abraham<sup>1</sup>, Joseph Mihaljevic<sup>1</sup>, Kris Murray<sup>2</sup>, and Christopher E. Doughty<sup>1</sup> 

<sup>1</sup> School of Informatics, Computing and Cyber Systems, Northern Arizona University, Flagstaff, USA.

<sup>2</sup> School of Public Health, Imperial College London, South Kensington Campus, London, UK.

The supplementary files needed for the successful running of the R/R markdown scripts (inc. previously published datasets) can be found at Zenodo (reserved DOI: 10.5281/zenodo.7521506).

Key figures:

----
<object data="https://github.com/Tomos/ExtinctHosts_PathogenRichness/blob/master/Figures/figure_process12.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://github.com/Tomos/ExtinctHosts_PathogenRichness/blob/master/Figures/figure_process12.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="http://yoursite.com/the.pdf">Download PDF</a>.</p>
    </embed>
</object>

![alt text](https://github.com/Tomos/ExtinctHosts_PathogenRichness/blob/master/Figures/figure_process12.pdf)

**Figure 1**: Schematic showing the generation of the extinction variables. All Phylacine range maps (left-hand side) represent extinct and extant species distributions under the assumption of no anthropogenic influence. Range maps were subdivided into the extant hosts (a), of which a subset was included in the Olival et al. (2017) study (b), as well as extinct species (c). Each of the host species from Olival et al. (2017) (b1) were compared to all other extant species (a1 – 4), to identify the number of extant sympatric species assuming no anthropogenic influence. For each sympatric species, the range map without anthropogenic pressure was compared to the corresponding map with anthropogenic pressure, to identify the number of local extinctions of extant species (represented by d1). The focal range map (b1) was also compared to all extinct ranges (c1 – 4), to identify all those that overlapped (c2 & c3). The Phylacine phylogeny was then used to determine the minimum, mean and median phylogenetic distances between the focal extant species (b1) and sympatric extinct species (c2 & c3). From all those that overlapped, the number of extinct species in the same family as the focal extant species were also determined. Animal silhouettes from PhyloPic.

---

![alt text](https://github.com/Tomos/ExtinctHosts_PathogenRichness/blob/master/Figures/global_local_extinction.png)

**Figure 2**: Maps generated from the Phylacine repository (version 1.2.1) showing the geographic richness of extinct mammal species (a) and the number of extant mammal species that are estimated to have been lost per pixel (local extinctions) due to anthropogenic influences (b), recorded since 130,000 ybp. Map (a) highlights hotspots of extinction in North and South America with no recorded mammal extinctions from central western Africa.

---

![alt text](https://github.com/Tomos/ExtinctHosts_PathogenRichness/blob/master/Figures/ShawGMPD2_GAMs-total-richness-bact-virus-helminth.png)

**Figure 3**: Partial dependence plots for GAMs fitted to viral (top), bacterial (middle) and helminth parasite (bottom) richness in extant mammalian species. Data on host-viral and -bacterial associations were taken from the Shaw et al. (2020) dataset, whereas the host-helminth data came from GMPD2. The GAMs were initially formulated by Olival et al. (2017), were subsequently used by Shaw et al. (2020), and included many of the same variables, including disease citations, host range and mammal sympatry. To these were added novel variables associated with Phylacine (version 1.2.1), summarized in Figure 1. Additionally, the first two mass-corrected principal components from Plourde et al. (2017) were also incorporated, which represent 85% of variation across 6 fast-slow life history traits.

---

