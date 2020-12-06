# 1. Set the Environment up ####

rm(list = ls(all = TRUE))
setwd("E:/PRIDE/R/SNA/Romania/Data")
#setwd("D:/PRIDE/R/SNA/Romania/Data")

getwd()

.libPaths('E:/PRIDE/R/Library')
#.libPaths('D:/PRIDE/R/Library')
.libPaths()

load('E:/PRIDE/R/SNA/Romania/Data/SNA.RData')
#load('D:/PRIDE/R/SNA/Romania/Data/SNA.RData')

# save(list = ls(all = TRUE), file='F:/PRIDE/R/SNA/Romania/Data/SNA.RData')
#save(list = ls(all = TRUE), file='D:/PRIDE/R/SNA/Romania/Data/SNA.RData')
save(list = ls(all = TRUE), file='E:/PRIDE/R/SNA/Romania/Data/SNA.RData')

library("sna")
library("igraph")
library("network")
library("GGally")
library("visNetwork")
library('extrafont')
library("reshape2")
library("tnet")
library("ggplot2")
library("plyr")
library("ggpubr")
library("Matrix")
library("DescTools")
library("nortest")
library("PMCMRplus")
library("sensitivity")
library("stringr")
library("splitstackshape")

citation(package="tnet")
citation(package= "igraph")
citation(package ="DescTools")
citation(package ="stats")
citation(package ="PMCMR")


# 2. Read the data ####

d <- read.csv(file = "Knowledge_transfer_SNA.Ro.csv", row.names = 1, header = T, as.is = T)  # First row contains row.names
d12 <- d 
colnames(d12) <- c(1:17)
rownames(d12) <- c(1:17)
nodes12 <- read.csv("nodes_all.csv", header=T, as.is=T)
nodes <- nodes12[, c(2,1,3,4)]

dm <- as.matrix(d)
dm_melt <- melt(dm)

colnames(dm_melt) <- c('snd','rec','weight')
dm_conf=subset(dm_melt,weight!=0)

a <- c(dm_conf$weight)
table(a) # 9 has been the case 12 times

(12*100)/75 # = 16%. So confirmation rate was 84 %. 

dm_melt_links=subset(dm_conf,weight!=9)
links <- dm_melt_links
write.csv(links, "links.all.directional.csv", row.names = F)

dm1 <- as.matrix(d12)
dm_melt1 <- melt(dm1)

colnames(dm_melt1) <- c('snd','rec','weight')
dm_melt_links1=subset(dm_melt1,weight!=0) #removes all rows of a certain value within a column
dm_melt_links1=subset(dm_melt_links1,weight!=9)
links12 <- dm_melt_links1
write.csv(links12, "links.all.directional12.csv", row.names = F)

links.1and2 <- links
links.1and2$weight <- ifelse(links.1and2$weight <= 2, 1,2)


#links <- read.csv("links.all.directional.csv", header = T, as.is = T)
#links12 <- read.csv("links.all.directional12.csv", header = T, as.is = T)

# Testing the preconditions for reconstruction ###

dm.test <- as.matrix(d12)
dm_melt.test <- melt(dm.test)
colnames(dm_melt.test) <- c('snd','rec','weight')
class(dm_melt.test)
dm_melt.test <- subset(dm_melt.test,snd!=10) #removes all the senders that are non-respondents 
dm_melt.test <- subset(dm_melt.test,snd!=11)
dm_melt.test <- subset(dm_melt.test,weight!=0)
dm_melt.test <- subset(dm_melt.test,weight!=9)

dm_nonres <- dm_melt.test[dm_melt.test$rec == 10 | dm_melt.test$rec == 11, ]
dm_res <- dm_melt.test[dm_melt.test$rec != 10 & dm_melt.test$rec != 11, ]

hist(dm_nonres$weight, breaks=10)
hist(dm_res$weight)
nonresCount <- as.matrix(table(dm_nonres$weight))
#nonresCount <- rbind(nonresCount, "4" = c(0))
nonresCount[,1] <- nonresCount[,1]/sum(nonresCount[,1])
resCount <- as.matrix(table(dm_res$weight))
resCount[,1] <- resCount[,1]/sum(resCount[,1])

# Shall we not have normalized the values by dividing them by number of respondents and Number of non-respondents respectively? 

GTest(x = nonresCount[,1],
      p = resCount[,1],
      correct='none') # p-value = 0.99... p>0.05 --> no difference in distributions 

GTest(x = nonresCount[c(2,1,3,4),1], # Order messed up! Just for testing p < 0.05, different distributions
      p = resCount[,1],
      correct='none')


# 3. Calculate Node statistics ####

# Degrees

DegW <- degree_w(links12, measure=c("degree","output"), type="out", alpha=1) # Weighted degree

head(DegW); dim(DegW)
#DegW <- DegW[, -1]
colnames(DegW) <- c( "id", "outdegree", "Woutdegree")
DegW.In <- degree_w(links12, measure=c("degree","output"), type = "in", alpha=1)
DegW.In
colnames(DegW.In) <- c("id" , "indegree", "Windegree")
head(DegW.In)

# Betweenness

btw <- betweenness_w(links12, directed = T, alpha = 1) # pretty representative
mean(btw) # 14.11
colnames(btw) <- c("id", "betweenness")
head(btw); dim(btw); str(btw)
colnames(btw)

node.stats <- data.frame(DegW, DegW.In, btw)
node.stats <- cbind(DegW, DegW.In[,2:3], btw[,2])
head(node.stats)
str(node.stats)
node.stats <- as.data.frame(node.stats)
colnames(node.stats)[6] <- "betweenness"
node.stats$Wdegree <- node.stats$Woutdegree + node.stats$Windegree
head(node.stats); dim(node.stats)
node.stats$Degree <- node.stats$outdegree + node.stats$indegree
head(node.stats); dim(node.stats)
node.stats <- node.stats[,c(1,4,2,6,5,3,8,7)]
node.stats$intieweight <- node.stats$Windegree/node.stats$indegree
node.stats$outtieweight <- node.stats$Woutdegree/node.stats$outdegree
write.csv (node.stats, "D:/PRIDE/PhD/Publication/SNA_Romania/Revisions 1/Results/Node.Statistics.directional.csv", row.names = F)
cor(node.stats[,2:3]) # Tie reciprocity = 0.38
# I was Here calculating the correlation coefficient withoug the Gov.

node.stats1 <- node.stats[-c(10,11,12,13,17), ]
cor(node.stats1[,2:3]) # Tie reciprocity = 0.79


quantile(node.stats$Degree) # 0%     25%     50%     75%     100% 
#  1      4       8       11      13   

quantile(node.stats$betweenness) # 0%         25%       50%       75%         100% 
# 0.00000  0.00000 15.00000    31.83333     87.33333 


# 4. Plot the data ####

net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
#ajmatx <- as_adjacency_matrix(net, attr="weight")
#ajmatx
#netaj <- graph.adjacency(ajmatx, mode = ("directed"), weighted = T, diag = FALSE)
#plot(netaj)


plot(net, edge.arrow.size=.4, edge.curved=.1)

# Generate colors based on Legal Status typI:
colrs <- c("deepskyblue4", "tomato", "gold1","darkgreen")
V(net)$color <- colrs[V(net)$Ls.type]
deg <- degree(net)
deg[1:17] <- node.stats$Wdegree
V(net)$size <- deg/3


E(net)$width <- E(net)$weight/2
plot(net, edge.color="orange", edge.arrow.size=.3, vertex.label.cex=1.2)

legend(x=-1.5, y=-0.9, c("Academic", "Governmental", "Non-governmental", "Protected Areas"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=3, cex=1, bty="n", ncol=1)


# Now we plot network with highlighted edges Strong and Weak


net12 <- graph_from_data_frame(d=links12, vertices=nodes12, directed=T)
net1.2 <- graph_from_data_frame(d=links.1and2, vertices=nodes, directed=T)
#ajmatx <- as_adjacency_matrix(net, attr="weight")
#ajmatx
#netaj <- graph.adjacency(ajmatx, mode = ("directed"), weighted = T, diag = FALSE)
#plot(netaj)

plot(net1.2)
plot(net1.2, edge.arrow.size=.4, edge.curved=.1)
plot(net1.2)

# Generate colors based on Legal Status typI:
colrs <- c("deepskyblue4", "tomato", "gold1","darkgreen")
V(net1.2)$color <- colrs[V(net1.2)$Ls.type]
deg <- degree(net)
deg[1:17] <- node.stats$Wdegree
V(net1.2)$size <- deg/2

colrsE <- c("gray60", "gray0") 
E(net1.2)$color <- colrsE[E(net1.2)$weight]
#E(net1.2)$width <- E(net1.2)$weight
plot(net1.2, edge.arrow.size=.2, vertex.label.cex=1.2)

legend(x=-1.5, y=-0.7, c("Academic", "Governmental", "Non-governmental", "Protected Areas"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=3, cex=1, bty="n", ncol=1)


# Now we plot network of only Strong relations

links.strong <- links.1and2
links.strong <- subset(links.strong,  weight == 2)

net.strong <- graph_from_data_frame(d=links.strong, vertices=nodes, directed=T)
#ajmatx <- as_adjacency_matrix(net, attr="weight")
#ajmatx
#netaj <- graph.adjacency(ajmatx, mode = ("directed"), weighted = T, diag = FALSE)
#plot(netaj)

plot(net.strong)
plot(net.strong, edge.arrow.size=.4, edge.curved=.1)

# Generate colors based on Legal Status typI:
colrs <- c("deepskyblue4", "tomato", "gold1","darkgreen")
V(net.strong)$color <- colrs[V(net.strong)$Ls.type]
deg <- degree(net)
deg[1:17] <- node.stats$Wdegree
V(net.strong)$size <- deg/3.2

colrsE <- c("gray60") 
E(net1.2)$color <- colrsE[E(net1.2)$weight]
#E(net1.2)$width <- E(net1.2)$weight
plot(net.strong, edge.arrow.size=.3, vertex.label.cex=1.2)

legend(x=-1.5, y=-0.9, c("Academic", "Governmental", "Non-governmental", "Protected Areas"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=3, cex=1, bty="n", ncol=1)


constraint(net.strong, nodes = V(net.strong), weights = NULL)

deg.str <- degree(net.strong)


dgstr <- c(deg.str)
dgstr <- as.data.frame(dgstr)
node.stats.str <- node.stats
node.stats.str$strong.links <- dgstr$dgstr

write.csv (node.stats.str, "D:/PRIDE/PhD/Publication/SNA_Romania/Revisions 1/Results/Node.Stats.str.csv", row.names = F)


######### Actual strong links not 1 and 2

links.strong1 <- links
links.strong1 <- subset(links.strong1,  weight > 2)

links.strong2 <- links12
links.strong2 <- subset(links.strong2, weight > 2)
dim(links.strong1); dim(links.strong2)


DegW.strong <- degree_w(links.strong2, measure=c("degree","output"), type="out", alpha=1) # Weighted degree

head(DegW.strong); dim(DegW.strong)
#DegW <- DegW[, -1]
colnames(DegW.strong) <- c( "id", "outdegree", "Woutdegree")
DegW.In.strong <- degree_w(links.strong2, measure=c("degree","output"), type = "in", alpha=1)
DegW.In.strong
colnames(DegW.In.strong) <- c("id" , "indegree", "Windegree")
head(DegW.In.strong)
degr.str <- cbind(DegW.strong, DegW.In.strong)
degr.str <- as.data.frame(degr.str)
degr.str$Degree <- degr.str$outdegree + degr.str$indegree
nodes.str <- nodes
nodes.str <- nodes.str[-c(15), ]

net.strong1 <- graph_from_data_frame(d=links.strong1, vertices=nodes.str, directed=T)
#ajmatx <- as_adjacency_matrix(net, attr="weight")
#ajmatx
#netaj <- graph.adjacency(ajmatx, mode = ("directed"), weighted = T, diag = FALSE)
#plot(netaj)

plot(net.strong1)
plot(net.strong1, edge.arrow.size=.4, edge.curved=.1)

# Generate colors based on Legal Status typI:
colrs <- c("deepskyblue4", "tomato", "gold1","darkgreen")
V(net.strong1)$color <- colrs[V(net.strong1)$Ls.type]
deg <- degree(net)
deg[1:17] <- node.stats$Wdegree
V(net.strong1)$size <- deg/3.2

colrsE <- c("gray0") 
E(net.strong1)$color <- colrsE[E(net.strong1)$weight]
#E(net1.2)$width <- E(net1.2)$weight
plot(net.strong1, edge.arrow.size=.3, vertex.label.cex=1.2)

legend(x=-1.5, y=-0.9, c("Academic", "Governmental", "Non-governmental", "Protected Areas"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=3, cex=1, bty="n", ncol=1)


Burtconstraint <- constraint(net.strong1, nodes = V(net.strong1), weights = NULL)
Burtconstraint <- as.matrix(Burtconstraint)
Burtconstraint <- as.data.frame(Burtconstraint)
write.csv(Burtconstraint, 'D:/PRIDE/PhD/Publication/SNA_Romania/Revisions 1/Results/Burt.csv')
#constraint(net1.2, nodes = V(net1.2), weights = NULL)

#Burtconstraint[-15, ] <- 1 # I am manually entering a very high constraint score
#Burtconstraint[-9, ] <- 1

burtweak <- constraint(net12, nodes = V(net12), weights = NULL)
burtweak <- as.matrix(burtweak)
burtweak <- as.data.frame(burtweak)
betweenness_w(links12, directed = T, alpha = 1)

Betweenness.strong <- betweenness_w(links.strong2, directed = T, alpha=1)
write.csv(Betweenness.strong, 'D:/PRIDE/PhD/Publication/SNA_Romania/Revisions 1/Results/Betweenness.strong.csv')
Burtconstraint1 <- Burtconstraint
Burtconstraint1$V2 <- Burtconstraint1$V1*100
Burtconstraint1$V2 <- round(Burtconstraint1$V2, digits = 0)

Betweenness.strongwithoutOC <- betweenness_w(links.strong1, directed=T, alpha=1)
dim(links.strong1); dim(links.strong2) # both 32, 3

quantile(Burtconstraint1$V2) # 0%     25%       50%         75%        100% 
# 25     32         43         74         100

quantile(Burtconstraint1$V1)

Betweenness.strong <- as.data.frame(Betweenness.strong)

quantile(Betweenness.strong$betweenness) #  0%       25%       50%       75%        100%
# 0.00000  0.00000    20        49          89 


# 5. Calculate Network Statistics ####

E(net) # - 63 links
V(net) # - 17 nodes
edge_density(net, loops = F) # 0.23: this means 23% of all the possible ties are present.
distance_w(links12, directed = TRUE)
mean(distance_w(links12, directed = TRUE), na.rm = T) # 2.2 This means that on average nodes are 2.2 steps, with average tie weight of 2.6, away from each other.
centr_degree(net, mode = ("all"), loops = F, normalized = T) # The graph level centrality index = 0.20
centr_degree(net, mode = ("in"), loops = F, normalized = T) # 0.34
centr_degree(net, mode = ("out"), loops = F, normalized = T) # 0.34
centr_betw(net, directed = TRUE, normalized = T) # for weak and strong ties included the centralization is 18 %
centr_betw(net.strong1, directed = TRUE, normalized = T) # For strong only ------ 21 %
mean(node.stats$Degree) # Mean degre = 7
node.stats$Tieweight <-  node.stats$Wdegree/node.stats$Degree
node.stats$Tieweight
mean(node.stats$Tieweight)

mean(links$weight) # 2.6
mean(node.stats$Degree) # 7


### Proportion of strong and weak relationships.

linksStrW <- links
head (linksStrW)
linksStrW$weight <- ifelse(linksStrW$weight <= 2, 1, 4)
dim(linksStrW)
sum(linksStrW$weight == 4) # 37
sum(linksStrW$weight == 1) # 26
(37*100/63) # ~ 59 % strong relationships
(26*100/63) # ~ 41 % Weak relationships

sum(links$weight == 4)
(9*100/65) # ~ 0nly 14 % very strog relationships


# 6. Count the Themes ####

Knowledge.Them <- read.csv(file = "Themes_Confirmed.Ro.csv", row.names = 1, header = T, as.is = T)
#mode(Knowledge.Them) <- "numeric"
#m <- as.matrix(Knowledge.Them)

#Them.category <- read.csv2(file = "Themes_Knowledge_transfer.Ro.csv", row.names = 1, header = T, as.is = T)

Them.Kn <- as.matrix(Knowledge.Them)
Them.list <- melt(Them.Kn)

colnames(Them.list) <- c('snd','rec','weight')
Themes =subset(Them.list,weight!=99)

Kn <- Knowledge.Them
colnames(Kn) <- c(1:17)
rownames(Kn) <- c(1:17)

Them.Kn1 <- as.matrix(Kn)
Them.list1 <- melt(Them.Kn1)

colnames(Them.list1) <- c('snd','rec','weight')
Themes1 =subset(Them.list1,weight!=99)

#######################################
dim(Themes); dim(Them.list) # 63   3; 289   3
Themes.check <- Themes
Themes.check$tie_strength <- links12$weight
themes.single <- Themes
#themes.single <- subset(themes.single, !weight %in% ";")
dim(themes.single)
#themes.single=paste(themes.single$themes, collapse = ";")
#themes.single <- str_split(themes.single$weight, ";")
thmn <- cSplit (themes.single, "weight", ";")

thmn$na_count <- apply(thmn, 1, function(x) sum(is.na(x)))
thmn$number <- (5-thmn$na_count)
links12
thmn$tie_strength <- links12$weight
them.sing.mult <- thmn
theme_N.and.strength <- table(thmn$number, thmn$tie_strength)
theme_N.and.strength <- t(theme_N.and.strength)

colnames(thmn) <-  c("snd","rec","them_id_1","them_id_2", "them_id_3", "them_id_4", "them_id_5", "bla","total_them_N", "tie_strength")
thmn <- thmn[ ,-8]
write.csv(thmn, "D:/PRIDE/PhD/Publication/SNA_Romania/Revisions 1/Results/Thems_and_Tie_strength.csv")


######### Test Theme and the tie strength association ######

shapiro.test(thmn$total_them_N) #p-value = 3.26e-09 - no normal ##chisq must be normal!
test<-wilcox.test(thmn$tie_strength,thmn$total_them_N)
test # Significantly different p-value = 0.0001504

########## A small test to see what happens if I exclude an outlier (link with 5 relational themes)
thmnTest <- thmn[-58, ]
thmnTest1 <- thmnTest [-55, ]
thmnTest2 <- thmnTest [-15, ]
thmnTest3 <- thmnTest2 [-15, ]
thmnTest4 <- thmnTest3 [-37, ]
thmnTest5 <- thmnTest4 [-52, ]
thmnTest6 <- thmnTest5 [-58, ]
shapiro.test(thmnTest$total_them_N) #p-value = 3.26e-09 - no normal ##chisq must be normal!
test.Test <- wilcox.test(thmnTest$tie_strength,thmnTest$total_them_N)
test.Test # Significantly different p-value 0.00009841

#
shapiro.test(thmnTest6$total_them_N) #p-value = 3.26e-09 - no normal ##chisq must be normal!
test.Test6 <- wilcox.test(thmnTest6$tie_strength,thmnTest6$total_them_N)
test.Test6 # Significantly different p-value = 0.0001418

#
shapiro.test(thmnTest4$total_them_N) #p-value = 3.26e-09 - no normal ##chisq must be normal!
test.Test4 <- wilcox.test(thmnTest4$tie_strength,thmnTest4$total_them_N)
test.Test4 # Significantly different p-value = 0.0001418


######### End of test

thmn$tie_strength <- as.numeric(thmn$tie_strength)
thmn$total_them_N <- as.numeric(thmn$total_them_N)

thmn[,8]
thmn$tie_strength <- ifelse(thmn$tie_strength >2,4,1)
m2 <- subset(thmn, tie_strength==4)
mean_2<-mean(m2$total_them_N)
m1 <- subset(thmn, tie_strength==1)
mean_1<-mean(m1$total_them_N)
boxplot(total_them_N~tie_strength,thmn)

data1 <- rbind(m1,m2)
class(data1)
data1 <- as.data.frame(data1)
colnames(data1) <-  c("snd","rec","them_id_1","them_id_2", "them_id_3", "them_id_4", "them_id_5", "No_contents_of_interaction", "Tie_strength")
data1$Tie_strength =factor(data1$Tie_strength)

means1 <- c(mean_1,mean_2)

ggplot (data1, aes(x =Tie_strength, y=No_contents_of_interaction)) + geom_boxplot() + scale_fill_grey() + theme_grey() +
  theme(text = element_text(size=20), axis.text.x = element_text(angle=0, hjust=1)) + 
  stat_summary(fun.y=mean, colour="black", geom="point", shape=18, size=3,show.legend = FALSE)


########## End of the test ###### Now below I will try to calculate individual stakeholder level themes of interaction

#Themes1 =subset(Them.list1,weight!=99)
#Themes.count <- Themes1
dim(links12) # 63   3
dim(Themes1) # 63   3
links.content <- links12
links.content.exp <- links.content[-1,]
nodes.exp <- nodes[-1,]
Themes.count.exp <- Themes.count[-1,]
#links12 <- links12[ ,-4]

links.content.exp$themes <- Themes.count.exp$weight
Them.categories <- read.csv(file = "Theme.categories.csv", header=T, as.is=T)

#links.content$exp <- paste(links.content$themes, collapse = ";")

i=8

table=data.frame()

table=c()


for (i in 2:17){
  x=subset(links.content.exp,snd==i)
  xn=nrow(x)
  y=subset(links.content.exp,rec==i)
  yn=nrow(y)
  name=nodes.exp$Organisation[nodes.exp$id==i]
  ls=nodes.exp$Legal.Status[nodes.exp$id==i]
  z=xn+yn
  degree=paste(z,"(",yn,",",xn,")",sep="")
  
  themey=paste(y$themes, collapse = ";")
  
  thy <- strsplit(themey, ";")[[1]]
  
  listy <- lapply(thy, function(x){Them.categories$category[Them.categories$id==x]})
  
  thy2 <- data.frame(t(count(unlist(listy, use.names=FALSE))))
  colnames(thy2) <- paste("in",thy2[1,])
  thy3 <- data.frame(thy2[-1,])
  row.names(thy3)=NULL
  
  themex=paste(x$themes, collapse = ";")
  
  thx <- strsplit(themex, ";")[[1]]
  
  listx <- lapply(thx, function(x){Them.categories$category[Them.categories$id==x]})
  
  thx2 <- data.frame(t(count(unlist(listx, use.names=FALSE))))
  class(thx2)
  colnames(thx2) <- paste("out",thx2[1,])
  thx3 <- data.frame(thx2[-1,])
  row.names(thx3)=NULL
  
  TAB=cbind(data.frame(name,ls,degree),thy3,thx3)
  table=rbind.fill(table,TAB)  
}
table[is.na(table)] <- 0

table2 <- table[,c(1,2,3,4,5,9,6,7,8)]
#write.csv(table2, "E:/PRIDE/PhD/Publication/SNA_Ukraine/Directed/Revisions_3/Manuscript/Results/themes.Org2.csv")
table3 <- melt(table2, id.var = c('name','ls','degree'))
deg <- c(node.stats$Degree)


############ my first looping experiment #######

j=8

table.exp=data.frame()

table.exp=c()


for (j in 2:17){
  x.snd=subset(links.content.exp,snd==j)
  x.snd.n=nrow(x.snd)
  x.rec=subset(links.content.exp,rec==j)
  x.rec.n=nrow(x.rec)
  name.x=nodes.exp$Organisation[nodes.exp$id==j]
  ls.x=nodes.exp$Legal.Status[nodes.exp$id==j]
  sum.x=x.snd.n+x.rec.n
  degree.x=paste(sum.x, "(",x.rec.n,",",x.snd.n,")",sep="")
  
  theme.x.rec=paste(x.rec$themes, collapse = ";")
  theme.x.rec.clean <- strsplit(theme.x.rec, ";")[[1]]
  theme.x.snd=paste(x.snd$themes, collapse = ";")
  theme.x.snd.clean <- strsplit(theme.x.snd, ";")[[1]]
  
  theme.x.rec.count <- as.matrix(table(theme.x.rec.clean))
  theme.x.rec.count <- t(theme.x.rec.count)
  colnames(theme.x.rec.count) <- paste("in", colnames(theme.x.rec.count))
  theme.x.snd.count <- as.matrix(table(theme.x.snd.clean))
  theme.x.snd.count <- t(theme.x.snd.count)
  colnames(theme.x.snd.count) <- paste("out", colnames(theme.x.snd.count))
  
  TAB.1=cbind(data.frame(name.x,ls.x,degree.x),theme.x.rec.count,theme.x.snd.count)
  table.exp=rbind.fill(table.exp,TAB.1)  
}
table.exp[is.na(table)] <- 0

table2.exp <- table.exp[,c(1,2,3,4,5,6,19,15,12,20,21,16,7,8,9,10,13,17,11,14,18)]
table2.exp.t <- t(table2.exp)
write.csv (table2.exp, "D:/PRIDE/PhD/Publication/SNA_Romania/Revisions 1/Results/Thems_and_Tie_strength.csv")


write.csv(table2.exp.t, "E:/PRIDE/PhD/Publication/SNA_Ukraine/Directed/Revisions_3/Manuscript/Results/themes.Org4.csv")
transp <- read.csv(file = 'themes.Org.transpose.csv', as.is = TRUE)
table3 <- melt(table2, id.var = c('name','ls','degree'))

################ Experiment ends here 


Themes.count <- Themes1
Th.all <- Themes.count$weight
Th.all <- paste(Th.all, collapse = ";")
Th.all1 <- strsplit(Th.all, ";")[[1]]
Th.all.clean <- gsub(" ", "", Th.all1)
table(Th.all.clean)
ab <- table(Th.all.clean)
class(ab)
ab1 <- as.data.frame(ab)
ab1[order(-ab1$Freq), ]

# Strong vs. Weak Themes 
Themes.count
links12
Themes.Str.vs.Weak <- Themes.count
Themes.Str.vs.Weak$Strength <- links12$weight
Themes.Str.vs.Weak$Strength <- ifelse(Themes.Str.vs.Weak$Strength <= 2, 1, 4)

sum(Themes.Str.vs.Weak$Strength == 4) # 37
sum(Themes.Str.vs.Weak$Strength == 1) # 26
Them.weak <- subset(Themes.Str.vs.Weak, Themes.Str.vs.Weak$Strength %in% 1)
Them.strong <- subset(Themes.Str.vs.Weak, Themes.Str.vs.Weak$Strength %in% 4)

th.weak <- Them.weak
th.weak <- th.weak[,-4]
th.weak <- th.weak$weight
th.weak <- paste(th.weak, collapse = ";")
th.weak1 <- strsplit(th.weak, ";")[[1]]
th.clean.weak <- gsub(" ", "", th.weak1)
table(th.clean.weak)
t.weak <- table(th.clean.weak)
class(t.weak)
t.weak1 <- as.data.frame(t.weak)
t.weak1[order(-t.weak1$Freq), ]


th.strong <- Them.strong$weight
th.strong <- paste(th.strong, collapse = ";")
th.strong1 <- strsplit(th.strong, ";")[[1]]
th.clean.strong <- gsub(" ", "", th.strong1)
table(th.clean.strong)
t.strong <- table(th.clean.strong)
class(t.strong)
t.strong1 <- as.data.frame(t.strong)
t.strong1[order(-t.strong1$Freq), ]

# SO from results above we see no associacion between the content of interactions and the strength of the tie.

######### Build a Pie Chart of Themes ########

# Theme Category Pie Chart
slices <- c(59,45)
lbls <- c("Collaboration relations", "Communication relations")
pie(slices, labels = lbls, col = terrain.colors(3), main="Pie Chart of identified themes")

legend("topleft", c("Collaboration relations", "Communication relations"), cex = 0.8,
       fill = terrain.colors(length(3)))


# Collaboration pie chart
slices.2 <- c(25, 13, 9, 7, 5)
lbls.2 <- c("Environmental projects", "Research", "Conservation planning", "Commercial fishing", "Sturgeon conservation")
pie(slices.2, labels = lbls.2, col = terrain.colors(6), main="Collaboration Relations")

# Communication pie chart
slices.1 <- c(18, 10, 10, 7)
lbls.1 <- c("Biodiversity data", "Environmental data", "Permit request", "Expert knowledge")
pie(slices.1, labels = lbls.1, col = terrain.colors(5), main="Communication Relations")


# 7. patterns according to Legal Status ####

# Governmental
head(nodes12)
nodes.Gov <- nodes12[nodes12$Ls.type == 2, ]
nodes.Gov
nodes.Gov.id <- c(nodes.Gov$id)
nodes.Gov.id

head(links12)
links.Gov.snd <- links12[links12$snd %in% nodes.Gov.id, ]
links.Gov.rec <- links12[links12$rec %in% nodes.Gov.id, ]
head(links.Gov.snd); dim(links.Gov.snd) # 15  3
head(links.Gov.rec); dim(links.Gov.rec) # 27  3

links.Gov.snd.rec <- links.Gov.snd[links.Gov.snd$rec %in% nodes.Gov.id, ]
links.Gov.snd.rec
colSums(links.Gov.snd.rec)[3] # Weight = 14
dim(links.Gov.snd.rec)[1] # 6

Governmental <- links.Gov.snd.rec
class(Governmental)

net.gov <- graph_from_data_frame(d=Governmental, vertices = nodes.Gov, directed = T)
plot(net.gov)
edge_density(net.gov, loops = F) # 0.3 

Gov12 <- Governmental
Gov12$weight <- ifelse(Gov12$weight <=2, 1,4)
Gov12Tab <- as.matrix(table(Gov12$weight)) # Within Governmental 4 weak and 2 strong relationships

# Academic

nodesAcad <- nodes12 [nodes12$Ls.type == 1, ]
Acadid <- c(nodesAcad$id)
linksAcadAcad <- subset(links12, links12$snd %in% Acadid & links12$rec %in% Acadid)
dim(linksAcadAcad)[1] # 19
net.Acad <- graph_from_data_frame(d=linksAcadAcad, vertices = nodesAcad, directed = T)
plot(net.Acad)
edge_density(net.Acad) # 0.26

Acad12 <- linksAcadAcad
Acad12$weight <- ifelse(Acad12$weight <=2, 1,4)
Acad12Tab <- as.matrix(table(Acad12$weight)) # Within academics 7 weak and 14 strong relationships

# NGO 

nodes.NGO <- nodes12[nodes12$Legal.Status == "Ngo", ]
nodes.NGO.id <- c(nodes.NGO$id)

links.NGO.snd <- links12[links12$snd %in% nodes.NGO.id, ]
links.NGO.rec <- links12[links12$rec %in% nodes.NGO.id, ]
head(links.NGO.snd); dim(links.NGO.snd) # 4  3
head(links.NGO.rec); dim(links.NGO.rec) # 5  3

links.NGO.snd.rec <- links.NGO.snd[links.NGO.snd$rec %in% nodes.NGO.id, ]
NGO <- links.NGO.snd.rec
net.NGO <- graph_from_data_frame(d=NGO, vertices = nodes.NGO, directed = T)
plot(net.NGO)
edge_density(net.NGO) # 0 

NGO12 <- NGO
NGO12$weight <- ifelse(NGO12$weight <=2, 1,4)
NGO12Tab <- as.matrix(table(NGO12$weight)) # zero relationships

# a. Statistics for different Sectors ####

degree_w(linksAcadAcad, measure=c("degree","output"), type="out", alpha=1) # Weighted degree
betweenness_w(linksAcadAcad, alpha = 1) # very interesting

mean(linksAcadAcad$weight) # 2.6
Governmental <- as.data.frame(Governmental)
mean(Governmental$weight) # 2.3
mean(NGO$weight) # NA
mean(PA$weight) # NA

# b. Network Statistics ####

V(net.Acad) # - 9 nodes
V(net.gov) # - 5 nodes
V(net.NGO) # - 3 nodes

E(net.Acad) # - 21 links
E(net.gov) # - 6 links 
E(net.NGO) # - 0 links


# Densities 

edge_density(net.Acad, loops = F) # 0.29 Acad - Acad
edge_density(net.gov, loops = F) # 0.3     Gov - Gov
edge_density(net.NGO, loops = F) # 0     NGO - NGO 


Acad.Gov.id <- c(1:13, 17)
nodesAcadGov <- nodes12 [nodes12$Ls.type == 1 | nodes12$Ls.type == 2 | nodes12$Ls.type == 4, ]
Acadid <- c(1:9)
Govid <- c(10:13, 17)

linksAcadGovsnd <- subset(links12, links12$snd %in% Acadid & links12$rec %in% Govid)
linksAcadGovrec <- subset(links12, links12$snd %in% Govid & links12$rec %in% Acadid)
dim(linksAcadGovsnd) # 19  3
dim(linksAcadGovrec) # 7  3
linksAcadGov <- rbind(linksAcadGovsnd, linksAcadGovrec)
dim(linksAcadGov)[1] # 26
net.AcadGov <- graph_from_data_frame(d=linksAcadGov, vertices = nodesAcadGov, directed = T)
plot(net.AcadGov)


Acad.GOv1 <- linksAcadGovsnd
Acad.GOv1$weight <- ifelse(Acad.GOv1$weight <=2, 1,4)
Acad.GOv1Tab <- as.matrix(table(Acad.GOv1$weight)) # 7 weak, 12 strong
Acad.GOv2 <- linksAcadGovrec
Acad.GOv2$weight <- ifelse(Acad.GOv2$weight <=2, 1,4)
Acad.GOv2Tab <- as.matrix(table(Acad.GOv2$weight)) # 4 weak, 3 strong

edge_density(net.AcadGov, loops = FALSE) # = 0.14
E(net.AcadGov) # 26 edges

# Acad Ngo

nodesAcadNGO <- nodes12 [nodes12$Ls.type == 1 | nodes12$Ls.type == 3, ]
NGOid <- c(14:16)
linksAcadNGOsnd <- subset(links12, links12$snd %in% Acadid & links12$rec %in% NGOid)
linksAcadNGOrec <- subset(links12, links12$snd %in% NGOid & links12$rec %in% Acadid)
dim(linksAcadNGOsnd) # 2 3
dim(linksAcadNGOrec) # 0 3 
linksAcadNGO <- rbind(linksAcadNGOsnd, linksAcadNGOrec)
net.AcadNGO <- graph_from_data_frame(d=linksAcadNGO, vertices = nodesAcadNGO, directed = T)
plot(net.AcadNGO)

edge_density(net.AcadNGO) # 0.015
10/(11*12) # = 0.07575758


nodesGovNGO <- nodes12 [nodes12$Ls.type == 2 | nodes12$Ls.type == 3, ]
linksGovNGOsnd <- subset(links12, links12$snd %in% Govid & links12$rec %in% NGOid)
linksGovNGOrec <- subset(links12, links12$snd %in% NGOid & links12$rec %in% Govid)
linksGovNGO <- rbind(linksGovNGOsnd, linksGovNGOrec)
net.GovNGO <- graph_from_data_frame(d=linksGovNGO, vertices = nodesGovNGO, directed = T)
plot(net.GovNGO)
edge_density(net.GovNGO) # 0.14

Gov.Ngo <- linksGovNGOsnd
Gov.Ngo$weight <- ifelse(Gov.Ngo$weight <=2, 1,4)
Gov.NgoTab <- as.matrix(table(Gov.Ngo$weight)) # 1 weak 3 strong

Ngo.Gov <- linksGovNGOrec
Ngo.Gov$weight <- ifelse(Ngo.Gov$weight <=2, 1,4)
Ngo.GovTab <- as.matrix(table(Ngo.Gov$weight)) # 1 weak 3 strong
E(net.GovNGO) # 8 edges

# Mean distance
mean(distance_w(linksAcadAcad), na.rm = T) # 2.3
mean(distance_w(Governmental, directed = T), na.rm = T) # 1.3
mean(distance_w(NGO, directed = T), na.rm = T) # NaN
mean(distance_w(PA, directed = F), na.rm = T) # NaN

DA <- degree_w(linksAcadAcad, measure=c("degree","output"), type="out", alpha=1)
DAD <-  as.data.frame(DA)
mean(DAD$output) #  6.2
mean(DAD$degree) # 2.3

Governmental <- as.matrix(Governmental)
mean(degree(net.gov)) #2.4 
mean(degree(net.Acad)) #4.6
mean(degree(net.NGO)) #0

# 7. Homophily ####

# 1 adjacency matrix d and nodes

d12
nodes12
head(links12)
head(nodes12)
unique(nodes12$Legal.Status)

# Create empty df
a <- c("Acad Acad", "Gov Gov", "Ngo Ngo", "Acad Gov", "Acad NGO", "Gov NGO")
edge.density.random.df <- as.data.frame(setNames(replicate(6, numeric(0), simplify = F), a))
dim(edge.density.random.df)
head(edge.density.random.df)

# a. Begin random ####


i=1

for(i in 1:1000){
  
  nodes1 <- as.data.frame(sample(nrow(nodes12)))
  nodes1
  names(nodes1)[1] <- "id"
  head(nodes1)
  nodes1 <- cbind(nodes1, nodes12[,3:4])
  nodes1 <- nodes1[order(nodes1[,c('id')]),]
  nodes1
  
  # Governmental
  head(nodes1)
  nodes1.Gov <- nodes1[nodes1$Legal.Status == "Gov", ]
  nodes1.Gov
  nodes1.Gov.id <- c(nodes1.Gov$id)
  nodes1.Gov.id
  
  head(links12)
  links.Gov.snd.rec1 <- links12[links12$snd %in% nodes1.Gov.id & links12$rec %in% nodes1.Gov.id, ]
  links.Gov.snd.rec1
  colSums(links.Gov.snd.rec1)[3] # Weight = 34
  dim(links.Gov.snd.rec1)[1] # 12
  
  Governmental.random <- links.Gov.snd.rec1
  class(Governmental.random)
  
  net.gov.random <- graph_from_data_frame(d=Governmental.random, vertices = nodes1.Gov, directed = T)
  plot(net.gov.random)
  edge_density(net.gov.random) # 0.6
  
  head(edge.density.random.df)
  edge.density.random.df[i, c("Gov.Gov")] <- edge_density(net.gov.random)
  
}

#edge.density.random.df <- edge.density.random.df[ ,-7]
hist(edge.density.random.df[, c("Gov.Gov")])
sd (edge.density.random.df[, c("Gov.Gov")])
edge.density.random.df$Gov.Gov
sort(edge.density.random.df$Gov.Gov) # no significantly different

# actual Gov-Gov edge density is 0.3 

# Acad ##

i=1

for(i in 1:1000){
  
  nodes1 <- as.data.frame(sample(nrow(nodes12)))
  nodes1
  names(nodes1)[1] <- "id"
  head(nodes1)
  nodes1 <- cbind(nodes1, nodes12[,3:4])
  nodes1 <- nodes1[order(nodes1[,c('id')]),]
  
  
  # Academic
  head(nodes1)
  nodes1.Academic <- nodes1[nodes1$Legal.Status == "Acad", ]
  nodes1.Academic
  nodes1.Academic.id <- c(nodes1.Academic$id)
  nodes1.Academic.id
  
  head(links12)
  links.Academic.snd.rec1 <- links12[links12$snd %in% nodes1.Academic.id & links12$rec %in% nodes1.Academic.id, ]
  links.Academic.snd.rec1
  colSums(links.Academic.snd.rec1)[3] # Weight = 20
  dim(links.Academic.snd.rec1)[1] # 8
  
  Academic.random <- links.Academic.snd.rec1
  class(Academic.random)
  
  net.Academic.random <- graph_from_data_frame(d=Academic.random, vertices = nodes1.Academic, directed = T)
  plot(net.Academic.random)
  edge_density(net.Academic.random) # 0.11
  
  head(edge.density.random.df)
  edge.density.random.df[i, c("Acad.Acad")] <- edge_density(net.Academic.random)
  
}

#edge.density.random.df <- edge.density.random.df[ ,-8]
head(edge.density.random.df)
hist(edge.density.random.df[, c("Acad.Acad")])
sd (edge.density.random.df[, c("Acad.Acad")])
sort(edge.density.random.df$'Acad.Acad')

# Actual Density is 0.29  

# NGO

i=1

for(i in 1:1000){
  
  nodes1 <- as.data.frame(sample(nrow(nodes12)))
  nodes1
  names(nodes1)[1] <- "id"
  head(nodes1)
  nodes1 <- cbind(nodes1, nodes12[,3:4])
  nodes1 <- nodes1[order(nodes1[,c('id')]),]
  
  
  # NGO
  head(nodes1)
  nodes1.NGO <- nodes1[nodes1$Legal.Status == "Ngo", ]
  nodes1.NGO
  nodes1.NGO.id <- c(nodes1.NGO$id)
  nodes1.NGO.id
  
  head(links12)
  links.NGO.snd.rec1 <- links12[links12$snd %in% nodes1.NGO.id & links12$rec %in% nodes1.NGO.id, ]
  links.NGO.snd.rec1
  colSums(links.NGO.snd.rec1)[3] # Weight = 3
  dim(links.NGO.snd.rec1)[1] # 14
  
  NGO.random <- links.NGO.snd.rec1
  class(NGO.random)
  
  net.NGO.random <- graph_from_data_frame(d=NGO.random, vertices = nodes1.NGO, directed = T)
  plot(net.NGO.random)
  edge_density(net.NGO.random) # 0
  
  head(edge.density.random.df)
  edge.density.random.df[i, c("Ngo.Ngo")] <- edge_density(net.NGO.random)
  
}

head(edge.density.random.df)
hist(edge.density.random.df[, c("Ngo.Ngo")])
hist(edge.density.random.df[, c("Ngo.Ngo")], breaks=c(unique(edge.density.random.df$Ngo.Ngo)))
sd (edge.density.random.df[, c("Ngo.Ngo")])

# Density = 0
sort(edge.density.random.df$Ngo.Ngo) # significantly low
sort(unique(edge.density.random.df$Ngo.Ngo))


# b. Between Groups random ####

# Create empty df

i=1

for(i in 1:1000){
  
  nodes2 <- as.data.frame(sample(nrow(nodes12)))
  nodes2
  names(nodes2)[1] <- "id"
  head(nodes2)
  nodes2 <- cbind(nodes2, nodes12[,3:4])
  nodes2 <- nodes2[order(nodes2[,c('id')]),]
  nodes2
  
  # Acad - Gov
  head(nodes2)
  nodes2.Acad <- nodes2[nodes2$Legal.Status == "Acad", ]
  nodes2.Gov <- nodes2[nodes2$Legal.Status == "Gov", ]
  nodes2.Acad.id <- c(nodes2.Acad$id)
  nodes2.Gov.id <- c(nodes2.Gov$id)
  nodesAcadGov.random <- rbind(nodes2.Acad, nodes2.Gov)
  
  head(links12)
  links.AcadGov.snd <- subset(links12, links12$snd %in% nodes2.Acad.id & links12$rec %in% nodes2.Gov.id)
  links.AcadGov.rec <- subset(links12, links12$snd %in% nodes2.Gov.id & links12$rec %in% nodes2.Acad.id)
  dim(links.AcadGov.snd) # 5  3
  dim(links.AcadGov.rec) # 7 3
  links.AcadGov.random <- rbind(links.AcadGov.snd, links.AcadGov.rec)
  AcadGov.random <- links.AcadGov.random
  class(AcadGov.random)
  
  net.AcadGov.random <- graph_from_data_frame(d=AcadGov.random, vertices = nodesAcadGov.random, directed = T)
  plot(net.AcadGov.random)
  dim(links.AcadGov.random)
  
  edge_density(net.AcadGov.random) 
  
  head(edge.density.random.df)
  edge.density.random.df[i, c("Acad.Gov")] <- edge_density(net.AcadGov.random)
  
}


hist(edge.density.random.df[, c("Acad.Gov")])
sd (edge.density.random.df[, c("Acad.Gov")])
edge.density.random.df$Acad.Gov
# Actual density is 0.14
sort(edge.density.random.df$Acad.Gov)

#Acad-NGO

i=1

for(i in 1:1000){
  
  nodes3 <- as.data.frame(sample(nrow(nodes12)))
  nodes3
  names(nodes3)[1] <- "id"
  head(nodes3)
  nodes3 <- cbind(nodes3, nodes12[,3:4])
  nodes3 <- nodes3[order(nodes3[,c('id')]),]
  nodes3
  
  # Acad - NGO
  head(nodes3)
  nodes3.Acad <- nodes3[nodes3$Legal.Status == "Acad", ]
  nodes3.NGO <- nodes3[nodes3$Legal.Status == "Ngo", ]
  nodes3.Acad.id <- c(nodes3.Acad$id)
  nodes3.NGO.id <- c(nodes3.NGO$id)
  nodesAcadNGO.random <- rbind(nodes3.Acad, nodes3.NGO)
  
  head(links12)
  links.AcadNGO.snd <- subset(links12, links12$snd %in% nodes3.Acad.id & links12$rec %in% nodes3.NGO.id)
  links.AcadNGO.rec <- subset(links12, links12$snd %in% nodes3.NGO.id & links12$rec %in% nodes3.Acad.id)
  links.AcadNGO.random <- rbind(links.AcadNGO.snd, links.AcadNGO.rec)  
  colSums(links.AcadNGO.random)[3] # Weight = 38
  dim(links.AcadNGO.random)[1] # 15
  
  AcadNGO.random <- links.AcadNGO.random
  class(AcadNGO.random)
  
  net.AcadNGO.random <- graph_from_data_frame(d=AcadNGO.random, vertices = nodesAcadNGO.random, directed = T)
  plot(net.AcadNGO.random)
  
  densAcadNGOrandom <- edge_density(net.AcadNGO.random); densAcadNGOrandom
  
  head(edge.density.random.df)
  edge.density.random.df[i, c("Acad.NGO")] <- densAcadNGOrandom
  
}

hist(edge.density.random.df[, c("Acad.NGO")])
sd (edge.density.random.df[, c("Acad.NGO")])
edge.density.random.df$Acad.NGO
# edge density is 1.5%
sort(edge.density.random.df$Acad.NGO) # Significantly different from random expectation


# Gov - NGO

i=1

for(i in 1:1000){
  
  nodes3 <- as.data.frame(sample(nrow(nodes12)))
  nodes3
  names(nodes3)[1] <- "id"
  head(nodes3)
  nodes3 <- cbind(nodes3, nodes12[,3:4])
  nodes3 <- nodes3[order(nodes3[,c('id')]),]
  nodes3
  
  head(nodes3)
  nodes3.Gov <- nodes3[nodes3$Legal.Status == "Gov", ]
  nodes3.NGO <- nodes3[nodes3$Legal.Status == "Ngo", ]
  nodes3.Gov.id <- c(nodes3.Gov$id)
  nodes3.NGO.id <- c(nodes3.NGO$id)
  nodesGovNGO.random <- rbind(nodes3.Gov, nodes3.NGO)
  
  head(links12)
  links.GovNGO.snd <- subset(links12, links12$snd %in% nodes3.Gov.id & links12$rec %in% nodes3.NGO.id)
  links.GovNGO.rec <- subset(links12, links12$snd %in% nodes3.NGO.id & links12$rec %in% nodes3.Gov.id)
  links.GovNGO.random <- rbind(links.GovNGO.snd, links.GovNGO.rec)
  colSums(links.GovNGO.random)[3] # Weight = 20
  dim(links.GovNGO.random)[1] # 7
  
  GovNGO.random <- links.GovNGO.random
  class(GovNGO.random)
  
  net.GovNGO.random <- graph_from_data_frame(d=GovNGO.random, vertices = nodesGovNGO.random, directed = T)
  plot(net.GovNGO.random)
  
  densGovNGOrandom <- edge_density(net.GovNGO.random); densGovNGOrandom
  
  head(edge.density.random.df)
  edge.density.random.df[i, c("Gov.NGO")] <- densGovNGOrandom
  
}

hist(edge.density.random.df[, c("Gov.NGO")])
sd (edge.density.random.df[, c("Gov.NGO")])
edge.density.random.df$Gov.NGO
# actual density is 14
sort(edge.density.random.df$Gov.NGO)


# Test normality

ggdensity(edge.density.random.df[, c("Gov")], xlab = "Network Density")
ggqqplot(edge.density.random.df[, c("Gov")])
shapiro.test(edge.density.random.df[, c("Gov")]) # p-value = 0.0001238. the test is significant so distribution is not normal. 


# Checking Themes according to Sectors ####

#Academic - Academic

nodesAcad <- nodes12 [nodes12$Ls.type == 1, ]
Acadid <- c(nodesAcad$id)
Themes.AcadAcad <- subset(Themes1, Themes1$snd %in% Acadid & Themes1$rec %in% Acadid)
th1 <-Themes.AcadAcad$weight
th1 <- paste(th1,collapse = ";")
s <- strsplit(th1, ";" )[[1]]
s.clean <- gsub(" ", "", s)
table(s.clean)
a <- table(s.clean)
class(a)
a1 <- as.data.frame(a)
a1[order(-a1$Freq), ]

# Governmental - Governmental

nodesGov <- nodes12 [nodes12$Ls.type == 2, ]
Govid <- c(nodesGov$id)
Themes.GovGov <- subset(Themes1, Themes1$snd %in% Govid & Themes1$rec %in% Govid)
th2 <-Themes.GovGov$weight
th2 <- paste(th2,collapse = ";")
G <- strsplit(th2, ";" )[[1]]
G.clean <- gsub(" ", "", G)
table(G.clean)
Govv <- table(G.clean)
class(Govv)
Gov1 <- as.data.frame(Govv)
Gov1[order(-Gov1$Freq), ]

# Check themes Between stakeholder groups ####

nodesNGO <- nodes12 [nodes12$Ls.type == 3, ]
NGOid <- c(nodesNGO$id)

nodesPA <- nodes12 [nodes12$Ls.type == 4, ]
PAid <- c(nodesPA$id)

# Inter-sectoral Themes ####
# Acad - Gov
ThemesAcadGovsnd <- subset(Themes1, Themes1$snd %in% Acadid & Themes1$rec %in% Govid)
ThemesAcadGovrec <- subset(Themes1, Themes1$snd %in% Govid & Themes1$rec %in% Acadid)
dim(ThemesAcadGovsnd) # 19  3
dim(ThemesAcadGovrec) # 7  3
ThemesAcadGov <- rbind(ThemesAcadGovsnd, ThemesAcadGovrec)
th5 <-ThemesAcadGov$weight
th5 <- paste(th5,collapse = ";")
AcadG <- strsplit(th5, ";" )[[1]]
Acad.G.clean <- gsub(" ", "", AcadG)
table(Acad.G.clean)
AGv <- table(Acad.G.clean)
class(AGv)
AG1 <- as.data.frame(AGv)
AG1[order(-AG1$Freq), ]

sum(AG1$Freq) # = 60

#Acad - Ngo
ThemesAcadNGOsnd <- subset(Themes1, Themes1$snd %in% Acadid & Themes1$rec %in% NGOid)
ThemesAcadNGOrec <- subset(Themes1, Themes1$snd %in% NGOid & Themes1$rec %in% Acadid)
dim(ThemesAcadNGOsnd) # 2  3
dim(ThemesAcadNGOrec) # 0  3
ThemesAcadNGO <- rbind(ThemesAcadNGOsnd, ThemesAcadNGOrec)
th6 <-ThemesAcadNGO$weight
th6 <- paste(th6,collapse = ";")
AcadNg <- strsplit(th6, ";" )[[1]]
Acad.Ng.clean <- gsub(" ", "", AcadNg)
table(Acad.Ng.clean)
ANgo <- table(Acad.Ng.clean)
class(ANgo)
ANgo1 <- as.data.frame(ANgo)
ANgo1[order(-ANgo1$Freq), ]



# Gov- NGO
ThemesGovNGOsnd <- subset(Themes1, Themes1$snd %in% Govid & Themes1$rec %in% NGOid)
ThemesGovNGOrec <- subset(Themes1, Themes1$snd %in% NGOid & Themes1$rec %in% Govid)
dim(ThemesGovNGOsnd) # 6  3
dim(ThemesGovNGOrec) # 6  3
ThemesGovNGO <- rbind(ThemesGovNGOsnd, ThemesGovNGOrec)
th9 <-ThemesGovNGO$weight
th9 <- paste(th9,collapse = ";")
GovNGO <- strsplit(th9, ";" )[[1]]
Gov.NGO.clean <- gsub(" ", "", GovNGO)
table(Gov.NGO.clean)
GoNGO <- table(Gov.NGO.clean)
class(GoNGO)
GoNGO1 <- as.data.frame(GoNGO)
GoNGO1[order(-GoNGO1$Freq), ]


# 8. Sufficient - Yes or No ####

sufficient <-  read.csv(file = "Sufficiency.Ro.csv", row.names = 1, header = T)
suf <- as.matrix(sufficient)
suf1 <- melt(suf)
colnames(suf1) <- c("snd", "rec", "value")
suf1 <- subset(suf1, !suf1$value ==0)
dim(suf1) # 63 
sum(suf1$value == 1) # 24 sufficient relations
sum(suf1$value == 6) ## 20 Unknown

63-(24+20)# =19 insufficient
19+24 # =43 known sufficiency
(19*100/43) # = 44%

#hist(suf1$value)
suf.count <- as.matrix(table(suf1$value))
suf.count
sufc <- as.data.frame(suf.count)
colnames(sufc) <- c("Frequency")
#sufc$"Sufficiency" <- c("Sufficient", "Insufficient due to lack of funding", "Lack of interconnection", "Unknown", "Insufficient due to instability", "Political constraint", "Competitive organizations", "Information is difficult to exchange")
row.names(sufc) <- c(1:8)
sufc.fin <- sufc[c(1,2,8,6,3,5,7,4),]

sufficient12 <- sufficient
colnames(sufficient12) <- c(1:17)
rownames(sufficient12) <- c(1:17)
sufficient12

suf.values <- as.matrix(sufficient12)
suf.values1 <- melt(suf.values)
colnames(suf.values1) <- c("snd", "rec", "value")
suf.values1 <- subset(suf.values1, !suf.values1$value ==0)

sum(suf.values1$value == 1) # 24 sufficient relations
sum(suf.values1$value == 6) ## 20 Unknown


insuf <- table(suf.values1$value)
insuf <- as.data.frame(insuf)
insuf


#### Stakeholder categories and sufficiency of interactions ####
suf12 <- suf.values1
suf12 <- subset(suf.values1, !suf.values1$value ==1)
dim(suf12)
suf12 <- subset(suf12, !suf12$value ==6)

suf12
#links12$sndrec <- paste(links12$snd, links12$rec, sep = "")
#suf12$sndrec <- paste(suf12$snd, suf12$rec, sep = "")

nodesAcad
Acadid
sufAcadAcad <- subset(suf12, suf12$snd %in% Acadid & suf12$rec %in% Acadid)
dim(sufAcadAcad) # 7 3
Acadsuf <- as.matrix(table(sufAcadAcad$value))

Govid
sufGovGov <- subset(suf12, suf12$snd %in% Govid & suf12$rec %in% Govid)
dim(sufGovGov) # 0 3
Govsuf <- as.matrix(table(sufGovGov$value))

#Acad - Gov
sufAcadGovsnd <- subset(suf12, suf12$snd %in% Acadid & suf12$rec %in% Govid)
sufAcadGovrec <- subset(suf12, suf12$snd %in% Govid & suf12$rec %in% Acadid)
sufAcadGov <- rbind(sufAcadGovsnd, sufAcadGovrec)
dim(sufAcadGov) # 9  3
AcadGovsuf <- as.matrix(table(sufAcadGov$value))

# Acad - NGO

sufAcadNGOsnd <- subset(suf12, suf12$snd %in% Acadid & suf12$rec %in% NGOid)
sufAcadNGOrec <- subset(suf12, suf12$snd %in% NGOid & suf12$rec %in% Acadid)
sufAcadNGO <- rbind(sufAcadNGOsnd, sufAcadNGOrec)
dim(sufAcadNGO) # 1  3
AcadNGOsuf <- as.matrix(table(sufAcadNGO$value))

# Gov - NGO

sufGovNGOsnd <- subset(suf12, suf12$snd %in% Govid & suf12$rec %in% NGOid)
sufGovNGOrec <- subset(suf12, suf12$snd %in% NGOid & suf12$rec %in% Govid)
sufGovNGO <- rbind(sufGovNGOsnd, sufGovNGOrec)
dim(sufGovNGO) # 2  3
GovNGOsuf <- as.matrix(table(sufGovNGO$value))

# I stopped here - will return after I calculate sufficiency..  
#thmn.exp <-  thmn
#suf1
#thmn.exp <- subset(thmn.exp, sufficient!=1 & sufficient!=8)
#table(thmn.exp$sufficient, thmn.exp$total_them_N)

# which ties are insufficient?
thmn
suf1
thmn$insufficient <- suf1$value
thmn1 <- thmn
thmn1 <- subset(thmn1, !thmn1$insufficient ==1 & !thmn1$insufficient ==6)

Iamsleepy <- as.matrix(table(thmn1$insufficient, thmn1$tie_strength))

### count how many links have how many themes.
link.themN <- as.matrix(table(thmn$total_them_N))
link.themN1 <- as.matrix(table(thmn$total_them_N, thmn$insufficient))

# Now lets see whether the links with single theme are sufficient or not
thmn
suf1
dim(thmn); dim(suf1)
suf.known
thmn$snd %in% suf1$snd & thmn$rec %in% suf1$rec # well aligned - all TRUE.. so
thmn$sufficient <- suf1$value
table(thmn$total_them_N, thmn$sufficient)

