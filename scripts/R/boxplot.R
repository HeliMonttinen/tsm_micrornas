library(ggplot2)
library(RColorBrewer)

cols5=brewer.pal(5,name = "Set2")

dist = read.table("results/distance_unmasked_tsms_variants_line.txt")[,7]

dat = cbind.data.frame(set="tsm",dist)

for(i in 1:100){
    dist = read.table(paste0("results/genome_unmasked_corrected_tsms_distance/background_line_",i))[,7]
    dat = rbind.data.frame(dat,cbind.data.frame(set=paste0("bg_",i),dist))
}

dat$dist=abs(dat$dist)

bins <- c(seq(0,650,50))
dat2 <- c()
for(s in unique(dat$set)){
    tmp <- cut(subset(dat,set==s)$dist,breaks=c(-1,seq(0,250,50),500,1000,1500,2000,3000,4000,5000,1000000))
    tmp <- factor(tmp,levels = levels(tmp), labels = bins)
    dat2 <- rbind.data.frame(dat2,cbind.data.frame(set=s,counts=as.vector(table(tmp)),bins=bins))
}
print(dat2)
image <- ggplot()+
 geom_boxplot(data=subset(dat2,set!="tsm"),aes(x=bins,y=counts,group=bins),color=cols5[2],outlier.size=.5,size=0.25,outlier.shape=4)+
 geom_point(data=subset(dat2,set=="tsm"),aes(x=bins,y=counts),shape=21,size=2,fill=cols5[3])+
 ylim(-1,max(dat2$counts))+
 guides(color="none")+
 labs(x="Distance",y="Count")+
 scale_x_continuous(labels = c("0 bp","1-50 bp","51-100 bp","101-150 bp","151-200 bp","201-250 bp","<0.5 kb","<1 kb","<1.5 kb","<2 kb","<3 kb","<4 kb","< 5kb",">5 kb"),
                    breaks=c(seq(0,650,50)))+
 theme_classic()+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))

ggsave("results/transposons_enrichment/plot_line_tsm_variants.pdf", image, width=4.5,height=4.5)

#pdf("results/plot_line_tsm_variants.pdf")
#print(image)     # Plot 1 --> in the first page of PDF
#dev.off() 
