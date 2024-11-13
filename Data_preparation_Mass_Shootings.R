rm(list = ls())
graphics.off()
library(readxl)
library(multiUS)

Offender = read_excel("Mass_Shootings.xlsx", sheet = "Offender", skip = 1)
Crime = read_excel("Mass_Shootings.xlsx", sheet = "Crime", skip = 1)
Traumas = read_excel("Mass_Shootings.xlsx", sheet = "Traumas", skip = 1)
Crisis = read_excel("Mass_Shootings.xlsx", sheet = "Crisis", skip = 1)
Mental = read_excel("Mass_Shootings.xlsx", sheet = "Mental", skip = 1)
# Mental = Mental[,-9]
Social = read_excel("Mass_Shootings.xlsx", sheet = "Social", skip = 1)
Motivation = read_excel("Mass_Shootings.xlsx", sheet = "Motivation", skip = 1)
Victims = read_excel("Mass_Shootings.xlsx", sheet = "Victims", skip = 0, col_types = rep("numeric", 2))

Crime = rowMeans(Crime != "0")
Traumas = rowMeans(Traumas != "0")
Crisis = rowMeans(Crisis != "0")
Social = rowMeans(Social != "0")
Mental = rowMeans(Mental != "0")
Motivation = rowMeans(Motivation != "0")

VictimsAge = aggregate(Victims$Age, by = list(Victims$`Case #`), FUN = median, na.rm = T)
for(j in 1:max(Victims$`Case #`)) {
  if(any(is.na(Victims$Age[Victims$`Case #` == j]))) {
    Victims$Age[Victims$`Case #` == j][is.na(Victims$Age[Victims$`Case #` == j])] = VictimsAge[j,2]
  }
}
VictimsAge = VictimsAge[-c(146),2]

Offender$Race[Offender$Race != 0] = 1
Offender$`Military Service`[Offender$`Military Service` == 2] = 0
Offender$`Relationship Status`[Offender$`Relationship Status` == 2] = 1
Offender$`Relationship Status`[Offender$`Relationship Status` == 3] = 0

DATA = data.frame(Offender$`Number Killed`, Offender$`Number Injured`, Offender$`Total Firearms Brought to the Scene`, 
                  Offender$Age, VictimsAge, Offender$`Insider / Outsider`, Offender$Immigrant, Offender$`Relationship Status`,
                  Social, Crime, Traumas, Crisis, Mental, Motivation, row.names = Offender$`Case #`, check.names = T)

DATA = DATA[-145,] # Wilkinsburg shooting (March 9, 2016)
DATA = head(DATA, -5) # 1966-2022 (2023 is incomplete)
100 *  sum(is.na(DATA)) / prod(dim(DATA)) # % of missing data

names(DATA) = as.character(c("Killed", "Injured", "Firearms.brought.to.the.scene", 
                             "Age", "Victims.age", "Insider", "Immigrant", "Relationship.status",
                             "Social", "Crime", "Traumas", "Crisis", "Mental", "Motivation"))


DATA = KNNimp(data = DATA, k = round(sqrt(nrow(na.omit(DATA)))), meth = "median")
apply(DATA[,1:8], 2, table)
# DATA %>% mutate_at(vars(names(DATA)[1:8]), list(~ round(., 0)))
# DATA = DATA[DATA$.Number.Injured. < 866,]

n = nrow(DATA)
p = ncol(DATA)
pc = ncol(DATA)


# Summary statistics - Tables ---------------------------------------------

summary.stats = function(data) {
  tau = c(0.25, 0.50, 0.75)
  
  data = as.matrix(data)
  n = nrow(data)
  pc = ncol(data)
  
  nonbinary = (1:pc)[which(apply(data, 2, function(x) length(unique(x))) != 2)]
  data_nonbinary = data[,nonbinary]
  data_binary = data[,-nonbinary]
  
  out_nonbinary = round(t(apply(data_nonbinary, 2, function(x) c(min(x), quantile(x, probs = tau)[1], mean(x), quantile(x, probs = tau)[2], 
                                                                 quantile(x, probs = tau)[3], max(x)))), 3)
  colnames(out_nonbinary) = c("Minimum", "First quartile", "Mean", "Median", "Third quartile", "Maximum")
  
  out_binary = round(t(apply(data_binary, 2, function(x) c(sum(x == 1), 100 * mean(x)))), 3)
  colnames(out_binary) = c("Frequency", "Proportion (%)")
  
  return(list(out_nonbinary = out_nonbinary, out_binary = out_binary))
}
out.stats = summary.stats(DATA)
out.stats

library(xtable)
xtable(out.stats$out_nonbinary, digits = 2)
xtable(out.stats$out_binary, digits = 2)

# Plots
Year = head(Offender$Year, -5)
Year = Year[-145]
Deaths  = aggregate(DATA$Killed, by = list(Year), FUN = sum) # Fatalities per year
Injured  = aggregate(DATA$Injured, by = list(Year), FUN = sum) # Injured per year
No_shoot = aggregate(Year, by = list(Year), FUN = length) # Number of shootings per year

par(mai = c(1, 1, 0.3, 0.1), cex.lab = 2, cex.axis = 1.9, cex.main = 2)

library(zoo)
wnd = 5

plot(y = Deaths$x, x = Deaths$Group.1, type = "h", xlab = "Year", ylab = "Fatalities", lwd = 2)
lines(y = rollapply(Deaths$x, FUN = mean, width = wnd), x = tail(Deaths$Group.1, 54 - wnd + 1), type = "l", col = "red", lwd = 4)

# lines(y = No_shoot$x, x = Deaths$Group.1, type = "l", col = "blue", lwd = 2)
# lines(y = Deaths$x / No_shoot$x, x = Deaths$Group.1, type = "l", col = "red", lwd = 2)
# lines(y = Injured$x / No_shoot$x, x = Deaths$Group.1, type = "l", col = "blue", lwd = 2)

plot(y = No_shoot$x, x = Deaths$Group.1, type = "h", xlab = "Year", ylab = "Shootings", lwd = 2)
lines(y = rollapply(No_shoot$x, FUN = mean, width = wnd), x = tail(Deaths$Group.1, 54 - wnd + 1), type = "l", col = "red", lwd = 4)

