## ----------------------------------------------------------
##        PIQ footprint analysis for ATAC-seq peaks
## ----------------------------------------------------------

library(data.table)
library(ggplot2)
library(ggrepel)
library(viridis)
library(ggpubr)
library(latex2exp)

## We used PIQ to find footprints for 1,316 motifs for B6 ATAC-seq samples
## We summarized the FIMO outputs for B6 ATAC-seq samples and computed the aggregated Footprint Depth (FPD) for each motif.
load('data-analysis/ATAC-seq/mouse/RData/dt-PIQ-enrichment-analysis.RData')
dt <- dt[is.human | is.mouse] # 744

dt <- dt[nFp > 0 & nMotif > 0] # remove 2 motifs without footprints
dt$motif.name <- dt$label
dt$motif <- paste0(dt$motif, ';', dt$motif.name)
dt$label <- sapply(strsplit(dt$label, '_'), function(x)
  x[1])

#### compute some quantities
# mean(dt[motif.name %in% c('Pdx1', 'Mnx1', 'NFATC2')]$nFp.peak)
# mean(dt[type %in% c('beta', 'islet')][! motif.name %in% c('Pdx1', 'Mnx1', 'NFATC2', 'CTCF')]$nFp.peak)



# percentage of ATAC-seq peaks containing footprints
dt$PercAtac <- dt$count / 510.14

# pooling motifs for the same TF
dt2 <-
  dt[, .(type, mean(nFp), mean(Fp.purity), mean(PercAtac)), by = label]
dt2 <- dt2[!duplicated(dt2)]
colnames(dt2) <- c('label', 'type', 'nFp', 'Fp.purity', 'PercAtac')


motifs <- dt2[PercAtac >= 3][order(-PercAtac)]$label
motifs <- c(motifs, dt2[type != 'other']$label)

dt2[label == 'NFATC2']$type <- 'beta'
dt2[label == 'CTCF']$type <- 'other'

dt2$type <- as.character(dt2$type)
dt2[type == 'alpha']$type <- 'alpha-cell'
dt2[type == 'beta']$type <- 'beta-cell'
dt2[type == 'islet']$type <- 'alpha & beta'
dt2$type <-
  factor(dt2$type,
         levels = c('alpha-cell', 'beta-cell', 'alpha & beta', 'other'))

dt2$text.color <- 'black'
dt2[label %in% c('ZNF354C', 'Gata1', 'GATA3', 'Pax4')]$text.color <-
  'grey'
dt2$text.color <- factor(dt2$text.color, levels = c('black', 'grey'))

## % of ATAC peaks with Fp vs. nFp
dt2[!label %in% motifs]$label <- ''
dt2$size <- 1
dt2[label != '']$size <- 2
g <- ggplot(dt2, aes(x = log2(nFp + 1), y = PercAtac)) +
  geom_point(aes(
    alpha = 1,
    size = size,
    color = Fp.purity,
    shape = type
  )) +
  geom_label_repel(
    data = dt2,
    aes(label = label, fill = text.color),
    color = 'white',
    size = 3,
    # the font of point label
    force = 5,
    segment.color = "gray30",
    box.padding = unit(0.35, 'lines'),
    point.padding = unit(0.3, 'lines'),
    show.legend = F
  )
g <-
  g + labs(x = expression(paste(log[2], '(', 'Total number of footprints', ')')),
           y = 'Percentage of ATAC-seq peaks with footprints',
           color = expression(atop('Average TF', 'occupancy probability')))

g <- g + theme_classic() +
  scale_alpha(guide = 'none', range = c(0.9, 0.9)) +
  scale_color_viridis() +
  scale_size(guide = 'none', range = c(1, 3)) +
  scale_fill_manual(values = c('black' = 'black', 'grey' = 'grey')) +
  scale_shape_discrete(name = '',
                       labels = c(
                         expression(paste(alpha, '-cells')),
                         expression(paste(beta, '-cells')),
                         expression(paste(alpha, ' & ', beta)),
                         'other'
                       )) +
  theme(legend.position = c(0.25, 0.6),
        legend.title = element_blank())
#  annotate('text',x=4.7, y=9, size = 3, label='Average TF \n occupancy probability')


pdf('Figure-PIQ-enrichment.pdf',
    width = 12,
    height = 4)
g + theme(legend.position = 'none')
dev.off()

pdf('Figure-PIQ-enrichment-legend.pdf', onefile = F)
leg <- get_legend(g)
as_ggplot(leg)
dev.off()



# motifs <- c("MA0139.1;CTCF","PB0015.1;Foxa2_1",
#             "PB0041.1;Mafb_1","PB0042.1;Mafk_1","PB0119.1;Foxa2_2",
#             "PB0145.1;Mafb_2","PB0146.1;Mafk_2","PH0003.1;Arx",
#             "PH0039.1;Mnx1","PH0081.1;Pdx1","PH0082.1;Irx2",
#             "PH0111.1;Nkx2-2","PH0118.1;Nkx6-1_1","PH0119.1;Nkx6-1_2",
#             "PH0131.1;Pax4","PH0132.1;Pax6")
# motifs <- c(motifs, "MA0068.1;Pax4","MA0069.1;Pax6","MA0132.1;Pdx1","MA0047.2;Foxa2")

motifs <- c(
  "PB0042.1;Mafk_1",
  "PB0146.1;Mafk_2",
  "PH0131.1;Pax4",
  "PH0111.1;Nkx2-2",
  "PB0015.1;Foxa2_1",
  "PB0119.1;Foxa2_2",
  "PH0132.1;Pax6",
  "MA0132.1;Pdx1",
  "PH0118.1;Nkx6-1_1",
  "PH0119.1;Nkx6-1_2"
)

dt <- dt[order(-observed)]
all.index <-
  sapply(motifs, function(x) {
    grep(x, paste0(dt$motif, ';', dt$label))
  })

## list of similar TFs for each islet TF
sim.TFs <- vector('list', length(all.index))

## check the islet TF enrichment
CheckSimilarTF <- function(i,
                           width.offset = 1,
                           ic.offset = 0.1) {
  index <- all.index[i]
  width.ref <- dt$width[index]
  ic.ref <- dt$ic[index]
  res <-
    dt[abs(width - width.ref) <= width.offset][abs(ic - ic.ref) <= ic.offset]
  rank <- grep(dt$motif[index], res$motif)
  print(paste(dt$motif[index], dt$label[index], 'rank =', rank, 'out of', nrow(res)))
  return(res)
}

for (i in 1:length(all.index)) {
  sim.TFs[[i]] <-
    CheckSimilarTF(i, width.offset = 1, ic.offset = 0.2)
}



## Now perform the randomization test
max.iter <- 100000
obs.stat <- mean(dt[all.index]$observed)
# [1] 0.2959
set.seed(527)
count <- 0

randomDrawStatistic <- function() {
  ans <- 0
  for (i in 1:length(sim.TFs)) {
    dt.tmp <- sim.TFs[[i]]
    rand.index <- sample(1:nrow(dt.tmp), 1)
    ans <- ans + dt.tmp$observed[rand.index]
  }
  ans <- ans / length(sim.TFs)
  return(ans)
}

for (iter in 1:max.iter) {
  count <- count + as.integer(randomDrawStatistic() > obs.stat)
}

p.value <- count / max.iter
# [1] 0.01663
