load("./shared/all.de.txpts.deseq2.Rdata")
load("./shared/ddsTxi.deseq2.RData")

## Create fam indices ####
this.p <- .05
this.lfc <- 2
all.fams <- as.character(ddsTxi$fam_id)
u.fams <- unique(all.fams)
bio.state.rep <- as.character(ddsTxi$state)
fl.i <- which(bio.state.rep %in% "FL")
fl.i.fams <- unique(all.fams[fl.i])
nc.i <- which(bio.state.rep %in% "NC")
nc.i.fams <- unique(all.fams[nc.i])
sc.i <- which(bio.state.rep %in% "SC")
sc.i.fams <- unique(all.fams[sc.i])
ga.i <- which(bio.state.rep %in% "GA")
ga.i.fams <- unique(all.fams[ga.i])
wg.i <- which(bio.state.rep %in% "WG")
wg.i.fams <- unique(all.fams[wg.i])

# FL vs. GA ####
out.u <- combn(c(fl.i.fams,ga.i.fams),m = 2)
rm <- unique(c(which(out.u[1,] %in% fl.i.fams & out.u[2,] %in% fl.i.fams),
               which(out.u[1,] %in% ga.i.fams & out.u[2,] %in% ga.i.fams)))
out.u <- out.u[,-rm]
out.u <- unique(t(out.u))

fl.vs.ga <- lapply(1:nrow(out.u),function(z){
  the.fams <- c(out.u[z,1],out.u[z,2])
  the.first.fam.cols <- which(u.fams %in% the.fams[1])
  the.second.fam.cols <- which(u.fams %in% the.fams[2])
  all.comp <- c(the.first.fam.cols,the.second.fam.cols)
  this.min <- which.min(all.comp)
  f1 <- paste0(the.fams[this.min])
  f2 <- paste0(the.fams[-this.min])
  other.fam <- all.comp[-this.min] - all.comp[this.min]
  if(this.min == 2){
    fl.up <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange > this.lfc)))
    fl.down <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange < -this.lfc)))
  } else {
    fl.down <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange > this.lfc)))
    fl.up <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange < -this.lfc)))
    
  }
  list("up"=fl.up,
       "down"=fl.down)
})

# FL vs. NC ####
out.u <- combn(c(fl.i.fams,nc.i.fams),m = 2)
rm <- unique(c(which(out.u[1,] %in% fl.i.fams & out.u[2,] %in% fl.i.fams),
               which(out.u[1,] %in% nc.i.fams & out.u[2,] %in% nc.i.fams)))
out.u <- out.u[,-rm]
out.u <- unique(t(out.u))

fl.vs.nc <- lapply(1:nrow(out.u),function(z){
  the.fams <- c(out.u[z,1],out.u[z,2])
  the.first.fam.cols <- which(u.fams %in% the.fams[1])
  the.second.fam.cols <- which(u.fams %in% the.fams[2])
  all.comp <- c(the.first.fam.cols,the.second.fam.cols)
  this.min <- which.min(all.comp)
  f1 <- paste0(the.fams[this.min])
  f2 <- paste0(the.fams[-this.min])
  
  other.fam <- all.comp[-this.min] - all.comp[this.min]
  if(this.min == 2){
    fl.up <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange > this.lfc)))
    fl.down <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange < -this.lfc)))
  } else {
    fl.down <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange > this.lfc)))
    fl.up <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange < -this.lfc)))
    
  }
  list("up"=fl.up,
       "down"=fl.down)
})
# FL vs. SC ####
out.u <- combn(c(fl.i.fams,sc.i.fams),m = 2)
rm <- unique(c(which(out.u[1,] %in% fl.i.fams & out.u[2,] %in% fl.i.fams),
               which(out.u[1,] %in% sc.i.fams & out.u[2,] %in% sc.i.fams)))
out.u <- out.u[,-rm]
out.u <- unique(t(out.u))

fl.vs.sc <- lapply(1:nrow(out.u),function(z){
  the.fams <- c(out.u[z,1],out.u[z,2])
  the.first.fam.cols <- which(u.fams %in% the.fams[1])
  the.second.fam.cols <- which(u.fams %in% the.fams[2])
  all.comp <- c(the.first.fam.cols,the.second.fam.cols)
  this.min <- which.min(all.comp)
  f1 <- paste0(the.fams[this.min])
  f2 <- paste0(the.fams[-this.min])
  
  other.fam <- all.comp[-this.min] - all.comp[this.min]
  if(this.min == 2){
    fl.up <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange > this.lfc)))
    fl.down <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange < -this.lfc)))
  } else {
    fl.down <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange > this.lfc)))
    fl.up <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange < -this.lfc)))
    
  }
  list("up"=fl.up,
       "down"=fl.down)
})

# FL vs. WG ####
out.u <- combn(c(fl.i.fams,wg.i.fams),m = 2)
rm <- unique(c(which(out.u[1,] %in% fl.i.fams & out.u[2,] %in% fl.i.fams),
               which(out.u[1,] %in% wg.i.fams & out.u[2,] %in% wg.i.fams)))
out.u <- out.u[,-rm]
out.u <- unique(t(out.u))

fl.vs.wg <- lapply(1:nrow(out.u),function(z){
  the.fams <- c(out.u[z,1],out.u[z,2])
  the.first.fam.cols <- which(u.fams %in% the.fams[1])
  the.second.fam.cols <- which(u.fams %in% the.fams[2])
  all.comp <- c(the.first.fam.cols,the.second.fam.cols)
  this.min <- which.min(all.comp)
  f1 <- paste0(the.fams[this.min])
  f2 <- paste0(the.fams[-this.min])
  
  other.fam <- all.comp[-this.min] - all.comp[this.min]
  if(this.min == 2){
    fl.up <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange > this.lfc)))
    fl.down <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange < -this.lfc)))
  } else {
    fl.down <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange > this.lfc)))
    fl.up <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange < -this.lfc)))
    
  }
  list("up"=fl.up,
       "down"=fl.down)
})

# SC vs. NC ####
out.u <- combn(c(sc.i.fams,nc.i.fams),m = 2)
rm <- unique(c(which(out.u[1,] %in% sc.i.fams & out.u[2,] %in% sc.i.fams),
               which(out.u[1,] %in% nc.i.fams & out.u[2,] %in% nc.i.fams)))
out.u <- out.u[,-rm]
out.u <- unique(t(out.u))
sc.vs.nc <- lapply(1:nrow(out.u),function(z){
  the.fams <- c(out.u[z,1],out.u[z,2])
  the.first.fam.cols <- which(u.fams %in% the.fams[1])
  the.second.fam.cols <- which(u.fams %in% the.fams[2])
  all.comp <- c(the.first.fam.cols,the.second.fam.cols)
  this.min <- which.min(all.comp)
  f1 <- paste0(the.fams[this.min])
  f2 <- paste0(the.fams[-this.min])
  
  other.fam <- all.comp[-this.min] - all.comp[this.min]
  if(this.min == 2){
    fl.up <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange > this.lfc)))
    fl.down <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange < -this.lfc)))
  } else {
    fl.down <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange > this.lfc)))
    fl.up <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange < -this.lfc)))
    
  }
  list("up"=fl.up,
       "down"=fl.down)
})

# SC vs. GA ####
out.u <- combn(c(sc.i.fams,ga.i.fams),m = 2)
rm <- unique(c(which(out.u[1,] %in% sc.i.fams & out.u[2,] %in% sc.i.fams),
               which(out.u[1,] %in% ga.i.fams & out.u[2,] %in% ga.i.fams)))
out.u <- out.u[,-rm]
out.u <- unique(t(out.u))
sc.vs.ga <- lapply(1:nrow(out.u),function(z){
  the.fams <- c(out.u[z,1],out.u[z,2])
  the.first.fam.cols <- which(u.fams %in% the.fams[1])
  the.second.fam.cols <- which(u.fams %in% the.fams[2])
  all.comp <- c(the.first.fam.cols,the.second.fam.cols)
  this.min <- which.min(all.comp)
  f1 <- paste0(the.fams[this.min])
  f2 <- paste0(the.fams[-this.min])
  
  other.fam <- all.comp[-this.min] - all.comp[this.min]
  if(this.min == 2){
    fl.up <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange > this.lfc)))
    fl.down <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange < -this.lfc)))
  } else {
    fl.down <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange > this.lfc)))
    fl.up <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange < -this.lfc)))
    
  }
  list("up"=fl.up,
       "down"=fl.down)
})

# SC vs. WG ####
out.u <- combn(c(sc.i.fams,wg.i.fams),m = 2)
rm <- unique(c(which(out.u[1,] %in% sc.i.fams & out.u[2,] %in% sc.i.fams),
               which(out.u[1,] %in% wg.i.fams & out.u[2,] %in% wg.i.fams)))
out.u <- out.u[,-rm]
out.u <- unique(t(out.u))
sc.vs.wg <- lapply(1:nrow(out.u),function(z){
  the.fams <- c(out.u[z,1],out.u[z,2])
  the.first.fam.cols <- which(u.fams %in% the.fams[1])
  the.second.fam.cols <- which(u.fams %in% the.fams[2])
  all.comp <- c(the.first.fam.cols,the.second.fam.cols)
  this.min <- which.min(all.comp)
  f1 <- paste0(the.fams[this.min])
  f2 <- paste0(the.fams[-this.min])
  
  other.fam <- all.comp[-this.min] - all.comp[this.min]
  if(this.min == 2){
    fl.up <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange > this.lfc)))
    fl.down <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange < -this.lfc)))
  } else {
    fl.down <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange > this.lfc)))
    fl.up <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange < -this.lfc)))
    
  }
  list("up"=fl.up,
       "down"=fl.down)
})

# NC vs. WG ####
out.u <- combn(c(nc.i.fams,wg.i.fams),m = 2)
rm <- unique(c(which(out.u[1,] %in% nc.i.fams & out.u[2,] %in% nc.i.fams),
               which(out.u[1,] %in% wg.i.fams & out.u[2,] %in% wg.i.fams)))
out.u <- out.u[,-rm]
out.u <- unique(t(out.u))
nc.vs.wg <- lapply(1:nrow(out.u),function(z){
  the.fams <- c(out.u[z,1],out.u[z,2])
  the.first.fam.cols <- which(u.fams %in% the.fams[1])
  the.second.fam.cols <- which(u.fams %in% the.fams[2])
  all.comp <- c(the.first.fam.cols,the.second.fam.cols)
  this.min <- which.min(all.comp)
  f1 <- paste0(the.fams[this.min])
  f2 <- paste0(the.fams[-this.min])
  
  other.fam <- all.comp[-this.min] - all.comp[this.min]
  if(this.min == 2){
    fl.up <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange > this.lfc)))
    fl.down <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange < -this.lfc)))
  } else {
    fl.down <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange > this.lfc)))
    fl.up <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange < -this.lfc)))
    
  }
  list("up"=fl.up,
       "down"=fl.down)
})

# NC vs. GA ####
out.u <- combn(c(nc.i.fams,ga.i.fams),m = 2)
rm <- unique(c(which(out.u[1,] %in% nc.i.fams & out.u[2,] %in% nc.i.fams),
               which(out.u[1,] %in% ga.i.fams & out.u[2,] %in% ga.i.fams)))
out.u <- out.u[,-rm]
out.u <- unique(t(out.u))
nc.vs.ga <- lapply(1:nrow(out.u),function(z){
  the.fams <- c(out.u[z,1],out.u[z,2])
  the.first.fam.cols <- which(u.fams %in% the.fams[1])
  the.second.fam.cols <- which(u.fams %in% the.fams[2])
  all.comp <- c(the.first.fam.cols,the.second.fam.cols)
  this.min <- which.min(all.comp)
  f1 <- paste0(the.fams[this.min])
  f2 <- paste0(the.fams[-this.min])
  
  other.fam <- all.comp[-this.min] - all.comp[this.min]
  if(this.min == 2){
    fl.up <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange > this.lfc)))
    fl.down <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange < -this.lfc)))
  } else {
    fl.down <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange > this.lfc)))
    fl.up <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange < -this.lfc)))
    
  }
  list("up"=fl.up,
       "down"=fl.down)
})
# GA vs. WG ####
out.u <- combn(c(ga.i.fams,wg.i.fams),m = 2)
rm <- unique(c(which(out.u[1,] %in% ga.i.fams & out.u[2,] %in% ga.i.fams),
               which(out.u[1,] %in% wg.i.fams & out.u[2,] %in% wg.i.fams)))
out.u <- out.u[,-rm]
out.u <- unique(t(out.u))
ga.vs.wg <- lapply(1:nrow(out.u),function(z){
  the.fams <- c(out.u[z,1],out.u[z,2])
  the.first.fam.cols <- which(u.fams %in% the.fams[1])
  the.second.fam.cols <- which(u.fams %in% the.fams[2])
  all.comp <- c(the.first.fam.cols,the.second.fam.cols)
  this.min <- which.min(all.comp)
  f1 <- paste0(the.fams[this.min])
  f2 <- paste0(the.fams[-this.min])
  
  other.fam <- all.comp[-this.min] - all.comp[this.min]
  if(this.min == 2){
    fl.up <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange > this.lfc)))
    fl.down <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange < -this.lfc)))
  } else {
    fl.down <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange > this.lfc)))
    fl.up <- rownames(as.data.frame(subset(all.de.txpts[[all.comp[this.min]]][[other.fam]],padj < this.p & log2FoldChange < -this.lfc)))
    
  }
  list("up"=fl.up,
       "down"=fl.down)
})
######
all.res <- list(fl.vs.ga,fl.vs.nc,fl.vs.sc,fl.vs.wg,
                sc.vs.ga,sc.vs.nc,sc.vs.wg,
                nc.vs.ga,nc.vs.wg,
                ga.vs.wg)

these.txpts <- unique(unlist(lapply(fl.vs.ga,function(x) x$down)))
these.txpts <- unique(these.txpts,unique(unlist(lapply(fl.vs.ga,function(x) x$up))))
cts <- vst(ddsTxi[c(these.txpts),c(fl.i,ga.i)])
plotPCA(cts,intgroup="")

fl.vs.ga[[1]]$up
the.tab.up <- unique(unlist(lapply(fl.vs.ga,function(x)x$up)))
the.tab.down <- unique(unlist(lapply(fl.vs.ga,function(x)x$down)))

rm.fl <- unique(c(unlist(lapply(all.res[c(2:4)],function(x) 
  lapply(x,function(z) z$up)))))
rm.fl <- which(the.tab.up %in% rm.fl)
the.tab.up <- the.tab.up[-rm.fl]

rm.fl <- unique(c(unlist(lapply(all.res[c(2:4)],function(x) 
  lapply(x,function(z) z$down)))))
rm.fl <- which(the.tab.down %in% rm.fl)
the.tab.down <- the.tab.down[-rm.fl]

rm.ga <- unique(c(unlist(lapply(all.res[c(5,8)],function(x) 
  lapply(x,function(z) z$down))),unlist(lapply(all.res[c(10)],function(x) 
    lapply(x,function(z) z$up)))))
rm.ga <- which(the.tab.up %in% rm.ga)
the.tab.up <- the.tab.up[-rm.ga]

rm.ga <- unique(c(unlist(lapply(all.res[c(5,8)],function(x) 
  lapply(x,function(z) z$up))),unlist(lapply(all.res[c(10)],function(x) 
    lapply(x,function(z) z$down)))))
rm.ga <- which(the.tab.down %in% rm.ga)
the.tab.down <- the.tab.down[-rm.ga]

cts <- varianceStabilizingTransformation(ddsTxi[c(the.tab.up,the.tab.down),c(fl.i,ga.i)])
plotPCA(cts,intgroup="state")

### EXTRA ####
rm.ga <- unique(c(unlist(lapply(all.res[c(5,8)],function(x) 
  lapply(x,function(z) z$up))),unlist(lapply(all.res[c(10)],function(x) 
    lapply(x,function(z) z$down)))))

rm.ga <- which(the.tab.down %in% rm.ga)
the.tab2 <- the.tab.down[-rm.ga]

rm.fl <- unique(c(unlist(lapply(all.res[c(2:4)],function(x) 
  lapply(x,function(z) z$down)))))
rm.fl <- which(the.tab.down %in% rm.fl)
the.tab33 <- the.tab2[-rm.fl]

cts <- varianceStabilizingTransformation(ddsTxi[c(the.tab3,the.tab33),c(fl.i,ga.i)])
plotPCA(cts,intgroup="state")

all.res <- list(fl.vs.ga,fl.vs.nc,fl.vs.sc,fl.vs.wg,
                sc.vs.ga,sc.vs.nc,sc.vs.wg,
                nc.vs.ga,nc.vs.wg,
                ga.vs.wg)
unique(unlist(lapply(all.res[c(5,8)],function(x) 
  lapply(x,function(z) z$down))))

all.res.sub <- all.res[c(1:3,5,6,8)]
these.cts <- c()
for(each.comp in 1:length(all.res.sub)){
  all.unlist <- unique(unlist(lapply(all.res.sub[-each.comp],function(x) unlist(x))))
  rm.txts <- which(unlist(all.res.sub[each.comp]) %in% all.unlist)
  fam.rm <- unlist(all.res.sub[each.comp])[-rm.txts]
  the.tab <- sort(table(fam.rm),decreasing = T)
  these.cts <- unique(c(these.cts,names(the.tab)))
}

cts <- varianceStabilizingTransformation(ddsTxi[these.cts,which(ddsTxi$batch %in% "LGEP")])
plotPCA(cts,intgroup="state")


all.unlist <- unique(c(unlist(fl.vs.ga),unlist(fl.vs.nc),unlist(fl.vs.wg),unlist(sc.vs.nc),unlist(sc.vs.ga),unlist(sc.vs.wg)))
rm.txts <- which(unlist(fl.vs.sc) %in% all.unlist)
fl.sc.rm <- unlist(fl.vs.sc)[-rm.txts]
the.tab <- sort(table(fl.sc.rm),decreasing = T)
fl.sc.cts <- names(the.tab)

all.unlist <- unique(c(unlist(fl.vs.ga),unlist(fl.vs.sc),unlist(fl.vs.wg),unlist(sc.vs.nc),unlist(nc.vs.ga),unlist(nc.vs.wg)))
rm.txts <- which(unlist(fl.vs.nc) %in% all.unlist)
fl.nc.rm <- unlist(fl.vs.nc)[-rm.txts]
the.tab <- sort(table(fl.nc.rm),decreasing = T)
fl.nc.cts <- names(the.tab)

all.unlist <- unique(c(unlist(sc.vs.wg),unlist(fl.vs.sc),unlist(sc.vs.ga),unlist(fl.vs.nc),unlist(nc.vs.ga),unlist(nc.vs.wg)))
rm.txts <- which(unlist(sc.vs.nc) %in% all.unlist)
sc.nc.rm <- unlist(sc.vs.nc)[-rm.txts]
the.tab <- sort(table(sc.nc.rm),decreasing = T)
sc.nc.cts <- names(the.tab)


