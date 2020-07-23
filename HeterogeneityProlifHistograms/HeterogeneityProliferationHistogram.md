
Load packages. 

```{r}
require(ggplot2)
require(ggridges)
```

Load files and uniformize column names.

```{r}
uM5 <- read.table("AllFCS5uM.tsv", header = TRUE, sep = "\t", row.names = 1)
colnames(uM5) <- c("FITC", "Generation", "CellType")
uM5$log2FITC <- log2(uM5$FITC)
uM5$CombinedID <- paste(uM5$CellType, uM5$Generation)
uM5$Generation <- factor(uM5$Generation)

uM10 <- read.table("AllFCS10nM.tsv", header = TRUE, sep = "\t", row.names = 1)
colnames(uM10) <- c("FITC", "Generation", "CellType", "log2FITC", "CombinedID")
uM10$Generation <- factor(uM10$Generation)
```

Removal of -Inf values from the log2 column. 

```{r}
uM5 <- uM5[is.finite(uM5$log2FITC), ]
uM10 <- uM10[is.finite(uM10$log2FITC), ]
```

Plot generation. 

```{r}
plot5uM <- ggplot(uM5, aes(x = uM5$log2FITC, y = uM5$Generation, fill = uM5$CellType)) + geom_density_ridges(alpha = 0.7) + scale_fill_manual(values = c("#427028", "#7C285C"), name = "Cell type") + ylab("Number of doublings") + xlab("log2(FITC)") + ggtitle("Concentration: 5 uM") + theme_classic() + theme(axis.line = element_blank(), panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
plot5uM 

plot10uM <- ggplot(uM5, aes(x = uM5$log2FITC, y = uM5$Generation, fill = uM5$CellType)) + geom_density_ridges(alpha = 0.7) + scale_fill_manual(values = c("#427028", "#7C285C"), name = "Cell type") + ylab("Number of doublings") + xlab("log2(FITC)") + ggtitle("Concentration: 10 nM") + theme_classic() + theme(axis.line = element_blank(), panel.grid = element_blank(), plot.title = element_text(hjust = 0.5)) 
plot10uM 
```

Export and save plots. 

```{r}
ggsave(plot5uM, filename = "Distribution_5uM.pdf", width = 8, height = 7)
ggsave(plot10uM, filename = "Distribution_10uM.pdf", width = 8, height = 7)
```
