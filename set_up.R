#############################################################################
# SET-UP: CREATING ENCOUNTER HISTORIES FOR RUNNING IN MARK
#############################################################################

# Steps in this set-up file:
# A. Make the encounter data frame: Assigns site codes to the Filemaker 
#    encounter file and adds the land cover, outputting a merged data frame
#
# B. Subsetting the data
# 	1. Subset to only species included in analyses
# 	2. Subset to only sites with an encounter history of >1
# 	3. Subset to only individuals in which the sex is known
# 	4. Subset to AHY and ASY then remove the band_age column
# 	5. Remove outliers in mass and wing data
# 	6. Subset sites with at least 3 years of point count observations
#
# C. Add an "OBSERVATION TYPE" Column: This will identify sites as those 
#    in which participants report data(1) versus those in which they do not(0).
#
# D. Create an encounter history ("ch") column
#
# E. Create a body condition index column: Creates a linear model of mass~
#    wing chord + sex + year + julian day for each species. The residuals of 
#    this relationship(scaled from 0 to 1) is the body condition index.
#
# F. Determine inter- and intraspecific abundance as well as the relative
#    abundance per species
#
# G. Create data frames (by species) to be used for MARK analyses

#---------------------------------------------------------------------------
# Set the libraries to be used in this file:
#---------------------------------------------------------------------------

library(outliers)
library(reshape)


#============================================================================
# A. Make the encounter data frame
#============================================================================

# Read in files:

setwd('C:/Users/Brian/Documents/survivorship/mark_input_files')

e = read.csv('encounters.csv')
codes = read.csv('site_codes.csv')
lc = read.csv('pts.can.imp.6.6.csv')

setwd('C:/Users/Brian/Documents/survivorship')

# Combine banding data with site codes:

a = merge(e,codes,by.x =24,by.y=1, all.x = T, all.y = F)

# Combine banding data with land cover:

a = merge(a,lc,by.x=26,by.y=1, all.x = T, all.y = F)
a = a[!is.na(a$imp100),]

colnames(a)[1] = 'site'

# Remove unwanted columns:

a = a[,-c(2,19,21,22,27:37)]

# I want imp500 to be listed as proportional LC rather than %

a$imp500 = a$imp500/100


#============================================================================
# B. Subsetting the data
#============================================================================

#----------------------------------------------------------------------------
# 1. Remove species not analyzed
#----------------------------------------------------------------------------

# Convert species column to all upper-case:

a$species = toupper(a$species)

# Subset to only NN species (not including NOMO):

a = a[a$species == 'AMRO'|a$species == 'CACH'|a$species == 'CARW'|a$species == 'GRCA'|a$species == 'HOWR'|a$species == 'NOCA'|a$species == 'SOSP',]
a$species = factor(a$species)


#----------------------------------------------------------------------------
# 2. Remove sites that have never had a re-encounter
#----------------------------------------------------------------------------

# Determine the total # of encounters for a given individual:

b.hist = ifelse(a[4:16] == '.', 0,ifelse(a[4:16]==0,0,1))
b.hist = rowSums(b.hist)

b = data.frame(a, b.hist)

# Determine the maximum number of encounters per site:

d1 = data.frame(tapply(b$b.hist,b$site,max))
	d1 = data.frame(row.names(d1),d1[1])
	colnames(d1) = c('site','enc')

# Remove sites that have never had an encounter:

d1 = na.omit(d1)

# Remove sites that have only had 1 encounter:

d1 = d1[d1$enc>1,]

# Merge with the data frame, keeping only matching records:

a = merge(a, d1, all.x = F, all.y = T)

# Replace the max number of encounters per site with the number of encounters
# per individual:

b.hist  = ifelse(a[4:16] == '.', 0,ifelse(a[4:16]==0,0,1))

a$enc = rowSums(b.hist)

# There is at least one individual without a capture history (No Band) ... remove it:

a = a[a$enc>0,]

# Finally, make sure the factor levels in "a" match the site subset:

a$site = factor(a$site)

#----------------------------------------------------------------------------
# 3. Subset to only individuals in which the sex is known
#----------------------------------------------------------------------------

a$sex = as.factor(toupper(a$sex))
a = a[!is.na(a$sex),]
a = a[a$sex == 'M'|a$sex == 'F',]

#----------------------------------------------------------------------------
# 4. Subset to AHY and ASY then remove the band_age column
#----------------------------------------------------------------------------

a = a[a$band_age!='HY',]
a = a[,-3]

#----------------------------------------------------------------------------
# 5. Test for and remove outliers in mass and wing data
#----------------------------------------------------------------------------
# Note: Outliers were determined using Grubbs test to find limits

# 5a. Subset species function:

spec.sub = function(species){
	d = a[a$species == species,]
	d = d[!is.na(d$mass)&!is.na(d$wing),]
	d$mw = d$mass/d$wing
	d
	}

amro = spec.sub('AMRO')
cach = spec.sub('CACH')
carw = spec.sub('CARW')
grca = spec.sub('GRCA')
howr = spec.sub('HOWR')
noca = spec.sub('NOCA')
sosp = spec.sub('SOSP')


# 5b. Cycle through grubbs tests to identify and remove outlier(s):

# amro

amro = amro[amro$mw>.265&amro$mw<.809,]

grubbs.test(amro$mw)

hist(amro$mw)


# cach

cach = cach[cach$mw>0&cach$mw<.209,]

grubbs.test(cach$mw)

hist(cach$mw)

# carw

carw = carw[carw$mw>0&carw$mw<.46,]

grubbs.test(carw$mw)

hist(carw$mw)

# grca


grca = grca[grca$mw>.234&grca$mw<.54,]

grubbs.test(grca$mw)

hist(grca$mw)


# howr


howr = howr[howr$mw>.139&howr$mw<.29215,]

grubbs.test(howr$mw)

hist(howr$mw)


# noca

noca = noca[noca$mw>.26&noca$mw<.608,]

grubbs.test(noca$mw)

hist(noca$mw)


# sosp

sosp = sosp[sosp$mw>.2075&sosp$mw<.6,]

grubbs.test(sosp$mw)

hist(sosp$mw)


# 5c. Combine the files back into the primary data frame:

a = rbind(amro, cach, carw, grca, howr, noca, sosp)

# 5d. Remove the mass-to-wing ratio column

a = a[,-24]

#============================================================================
# C. Add an "OBSERVATION TYPE" Column
#============================================================================
# Note: this will identify sites as those in which participants report data(1)
# versus those in which they do not(0).

# Determine the number of observations by participants:

p.hist = ifelse(a[3:15]=='rp'|a[3:15]=='rprt'|a[3:15]=='rprc',1,0)
p.hist = rowSums(p.hist)

p1 = data.frame(a, p.hist)

# Determine the maximum number of participant encounters per site:

d1 = data.frame(tapply(p1$p.hist,p1$site,max))
	d1 = data.frame(row.names(d1),d1[1])
	colnames(d1) = c('site','participation')

# Remove sites that have never had an encounter:

d1 = na.omit(d1)

# Merge with the data frame, keeping only matching records:

a = merge(a, d1, all.x = F, all.y = T)

# Recode the participation column as active or inactive:

a$participation = factor(ifelse(a$participation>0, 'active', 'inactive'))

# Some of the marked birds had no band added ... remove them:

a = a[a$band_num != 'NB',]

#============================================================================
# D. Create an encounter history ("ch") column
#============================================================================

# Encounter history:
# "." = unsampled 
# "0" = unencountered
# "1" = encountered

ch1 = ifelse(a[3:15] == '.',0,1)

# As encounter history column:

a$ch = paste(ch1[,1],ch1[,2],ch1[,3],ch1[,4],ch1[,5],ch1[,6],ch1[,7],ch1[,8],ch1[,9],ch1[,10],ch1[,11],ch1[,12],ch1[,13],sep='')

# Remove encounter histories separated by years:

a = a[,-c(3:15)]


#============================================================================
# E. Create a body condition index column
#============================================================================

# A function to create a dataframe to analyze body condition:

bci.df = function(species){
	df1 = spec.sub(species)
	mass = df1$mass
	wing = df1$wing
	sex = df1$sex
	imp500 = df1$imp500
	date1 = strptime(df1$d_banded, "%m/%d/%Y")
	year = as.numeric(format(date1,format = '%Y'))
	day = as.numeric(format(date1,format = '%j'))
	df2 = data.frame(df1$band_num,mass, wing, sex, year, day, imp500)
	colnames(df2)[1] = 'band_num'
	df2
	}

# Create separate data frames by species:

amro = bci.df('AMRO')
cach = bci.df('CACH')
carw = bci.df('CARW')
grca = bci.df('GRCA')
howr = bci.df('HOWR')
noca = bci.df('NOCA')
sosp = bci.df('SOSP')

# Calculate the body condition index:
# Following Peig and Green 2009

bci.fun = function(species){
	df1 = species
	L0 = mean(df1$wing)
	Mi= df1$mass
	Li = df1$wing
	bsma = summary(lm(log(Mi)~log(Li)))$coefficients[2,1]/cor(log(Mi),log(Li))
	Mi*((L0/Li)^bsma)
	}

# Rescale each body condition index from 0 to 1:

amro$bci = rescaler(bci.fun(amro), type = 'range')
cach$bci = rescaler(bci.fun(cach), type = 'range')
carw$bci = rescaler(bci.fun(carw), type = 'range')
grca$bci = rescaler(bci.fun(grca), type = 'range')
howr$bci = rescaler(bci.fun(howr), type = 'range')
noca$bci = rescaler(bci.fun(noca), type = 'range')
sosp$bci = rescaler(bci.fun(sosp), type = 'range')

# Bind into a new data frame:

bci.df = rbind(amro, cach, carw, grca, howr, noca, sosp)

# Remove all columns except band number, bci and year:

bci.df = bci.df[,-c(2:4,6:7)]

# Merge with the primary data frame by band number:

a = merge(a, bci.df)

# Make sure no duplicates have been generated:

a = unique(a)


#============================================================================
# G. Code observations as "sampled" or "unsampled" (dot method) for each year
#============================================================================

a.frame = a[,c(1:2,11:12)]

string.list = strsplit(a.frame$ch,'')

df = data.frame(matrix(as.numeric(unlist(string.list)), ncol =13, byrow = T))
head(df)

af2 = data.frame(a.frame,df)

year.count.fun = function(year.column){
	y1 = data.frame(tapply(af2[[year.column]],af2$site,max))
	y1 = data.frame(row.names(y1),y1)
	colnames(y1) = c('site','sampled')
	y1
	}

y = year.count.fun('X1')
y$y01 = year.count.fun('X2')[,2]
y$y02 = year.count.fun('X3')[,2]
y$y03 = year.count.fun('X4')[,2]
y$y04 = year.count.fun('X5')[,2]
y$y05 = year.count.fun('X6')[,2]
y$y06 = year.count.fun('X7')[,2]
y$y07 = year.count.fun('X8')[,2]
y$y08 = year.count.fun('X9')[,2]
y$y09 = year.count.fun('X10')[,2]
y$y10 = year.count.fun('X11')[,2]
y$y11 = year.count.fun('X12')[,2]
y$y12 = year.count.fun('X13')[,2]
colnames(y)[2] = 'y00'


af2$y01 = y$y01[match(af2$site, y$site)]
af2$y02 = y$y02[match(af2$site, y$site)]
af2$y03 = y$y03[match(af2$site, y$site)]
af2$y04 = y$y04[match(af2$site, y$site)]
af2$y05 = y$y05[match(af2$site, y$site)]
af2$y06 = y$y06[match(af2$site, y$site)]
af2$y07 = y$y07[match(af2$site, y$site)]
af2$y08 = y$y08[match(af2$site, y$site)]
af2$y09 = y$y09[match(af2$site, y$site)]
af2$y10 = y$y10[match(af2$site, y$site)]
af2$y11 = y$y11[match(af2$site, y$site)]
af2$y12 = y$y12[match(af2$site, y$site)]

af3 = af2

af3$participation = ifelse(af3$participation == 'active',1,0)

af3[af3$site == 'MARRPETMD2',] 
#ifelse(af3[[y.year]] == 1,'0','.')

recode.fun = function(X.year,y.year){
	test = as.factor(ifelse(af3[[X.year]] == 1, '1',
	ifelse(af3$participation == 'active',0,
	ifelse(af3[[y.year]] == 0,'.',0))))
	test
	}

af3$y00 = af3[['X1']]
af3$y01 = recode.fun('X2','y01')
af3$y02 = recode.fun('X3','y02')
af3$y03 = recode.fun('X4','y03')
af3$y04 = recode.fun('X5','y04')
af3$y05 = recode.fun('X6','y05')
af3$y06 = recode.fun('X7','y06')
af3$y07 = recode.fun('X8','y07')
af3$y08 = recode.fun('X9','y08')
af3$y09 = recode.fun('X10','y09')
af3$y10 = recode.fun('X11','y10')
af3$y11 = recode.fun('X12','y11')
af3$y12 = recode.fun('X13','y12')


af3$ch2 = paste(af3$y00,af3$y01,af3$y02,af3$y03,af3$y04,af3$y05,af3$y06,af3$y07,af3$y08,af3$y09,af3$y10,af3$y11,af3$y12, sep = '')

a$ch = af3$ch2

#============================================================================
# H. Create data frames (by species) to be used for MARK analyses
#============================================================================

# Create the data frame:

a = data.frame(a[7:9],a[11:12],a[14])

# Add an impervious^2 column

a$imp2 = a$imp500^2

# Re-order the columns

a = a[,c(2,5,1,6,3,7,4)]

# Rename the columns:

colnames(a)[c(5,7)] = c('imp','part')

# Create the subsets by species:

spec.sub = function(species){
	d = a[a$species == species,]
	d[,-1]				# Removes the species name column
	}

amro = spec.sub('AMRO')
cach = spec.sub('CACH')
carw = spec.sub('CARW')
grca = spec.sub('GRCA')
howr = spec.sub('HOWR')
noca = spec.sub('NOCA')
sosp = spec.sub('SOSP')

###############################################################################
# 						STOP							
###############################################################################
