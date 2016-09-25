mark.height <- 0.01
v.sep <- 0.008
bar.width <- 0.45
  
test.fractions <- function(fractions) {
  marks.lines <- data.frame()
  fractions$pos <- fractions$Fraction * fractions$number
  
  for (l_sec in levels(fractions$sec_simple)) {
    s_1 <- subset(fractions, sec_simple == l_sec & is_domain == "domain")
    s_2 <- subset(fractions, sec_simple == l_sec & is_domain == "non-domain")
    print(s_1)
    print(s_2)
    res <- prop.test(c(s_1$pos[1], s_2$pos[1]), c(s_1$number[1], s_2$number[1]))
    print(res)
    
    if (res$p.value < 0.05) {
      stars <- "*"
      
      if (res$p.value < 0.01)
        stars <- "**"
      
      if (res$p.value < 0.001)
        stars <- "***"
      
      # Horizontal line
      y.min <- 
    }
  }
}

test.fractions(fractions.domains.buried)


test.table <- function(fractions) {
  res <- data.frame()
  fractions$pos <- fractions$Fraction * fractions$number
    
  for (domain_status in c("domain", "non-domain")) {
    for (l_1 in levels(fractions$sec_simple)) {
      for (l_2 in levels(fractions$sec_simple)) {
        if (l_1 <= l_2) {
          next
        }
        s_1 <- subset(fractions, sec_simple == l_1 & is_domain == domain_status)
        s_2 <- subset(fractions, sec_simple == l_2 & is_domain == domain_status)
        
        test_res <- prop.test(c(s_1$pos[1], s_2$pos[1]), c(s_1$number[1], s_2$number[1]))
        
        res <- rbind(res, data.frame(sec_simple_1=l_1, sec_simple_2=l_2, domain=domain_status, p=test_res$p.value))
      }
    }
  }
  
  res
}

test.table.buried <- test.table(fractions.domains.buried)
test.table.exposed <- test.table(fractions.domains.exposed)

test.table.buried$p.adj <- p.adjust(c(test.table.buried$p, test.table.exposed$p), method="BH", n=24)[1:12]
test.table.exposed$p.adj <- p.adjust(c(test.table.buried$p, test.table.exposed$p), method="BH", n=24)[13:24]

write.table(test.table.buried, file=paste(results, "test_table_domains_buried.tsv"), sep="\t", row.names=F)
write.table(test.table.exposed, file=paste(results, "test_table_domains_exposed.tsv"), sep="\t", row.names=F)
