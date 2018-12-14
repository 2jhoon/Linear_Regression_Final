data.path = 
  '/Users/jeonghoonlee/Desktop/SSU/4학년 2학기/회귀분석 2/기말프로젝트/Final_Project/dat_train.csv'
data = read.csv(data.path)
head(data)
str(data)

### 범주형 & 연속형 변수 구분

dis_val = c("Sex", "smoking", "HTN", "DM")
dis_data = data[dis_val]

con_val = c('age', 'TCHOL', 'hdl', 'LDL', 'hemoglobin', 'WBC',
            'neutrophil', 'lymphocyte')
con_data = data[con_val]

### 연속형 변수 탐색 시각화
### pair plot
pairs(con_data, lower.panel = panel.smooth, upper.panel = panel.cor)
par(mfrow = c(2,4))

### box plot
boxplot(con_data$age, xlab = 'age', boxwex = 1.5); boxplot(con_data$TCHOL, xlab = 'TCHOL', boxwex = 1.5)
boxplot(con_data$hdl, xlab = 'hdl', boxwex = 1.5); boxplot(con_data$LDL, xlab = 'LDL', boxwex = 1.5)
boxplot(con_data$hemoglobin, xlab = 'hemoglobin', boxwex = 1.5); boxplot(con_data$WBC, xlab = 'WBC', boxwex = 1.5)
boxplot(con_data$neutrophil, xlab = 'neutrophil', boxwex = 1.5); boxplot(con_data$lymphocyte, xlab = 'lymphocyte', boxwex = 1.5)
par(mfrow = c(1,1))

### 범주형 변수 탐색
table = table(dis_data)
par(mfrow = c(1,4))
plot(factor(dis_data$Sex, labels = c('Female', 'Male')), main = 'Sex')
plot(factor(dis_data$smoking, labels = c('False', 'True')), main = 'smoking')
plot(factor(dis_data$HTN, labels = c('False', 'True')), main = 'HTN')
plot(factor(dis_data$DM, labels = c('False', 'True')), main = 'DM')

CAD.group = data[data$CADGROUP==1, ]
plot(factor(CAD.group$Sex, labels = c('Female', 'Male')), main = 'CAD - True: Sex')
plot(factor(CAD.group$smoking, labels = c('False', 'True')), main = 'CAD - True: smoking')
plot(factor(CAD.group$HTN, labels = c('False', 'True')), main = 'CAD - True: HTN')
plot(factor(CAD.group$DM, labels = c('False', 'True')), main = 'CAD - True: DM')

par(mfrow = c(1,3))
CAD.Female = CAD.group[CAD.group$Sex==1, ]
plot(factor(CAD.Female$smoking, labels = c('False', 'True')), main = 'CAD - True: smoking')
plot(factor(CAD.Female$HTN, labels = c('False', 'True')), main = 'CAD - True: HTN')
plot(factor(CAD.Female$DM, labels = c('False', 'True')), main = 'CAD - True: DM')
sum(CAD.Female$smoking==1) / nrow(CAD.Female)
sum(CAD.Female$HTN==1) / nrow(CAD.Female)
sum(CAD.Female$DM==1) / nrow(CAD.Female)

CAD.Male = CAD.group[CAD.group$Sex==2, ]
plot(factor(CAD.Male$smoking, labels = c('False', 'True')), main = 'CAD - True: smoking')
plot(factor(CAD.Male$HTN, labels = c('False', 'True')), main = 'CAD - True: HTN')
plot(factor(CAD.Male$DM, labels = c('False', 'True')), main = 'CAD - True: DM')
sum(CAD.Male$smoking==1) / nrow(CAD.Male)
sum(CAD.Male$HTN==1) / nrow(CAD.Male)
sum(CAD.Male$DM==1) / nrow(CAD.Male)

par(mfrow = c(1,2))
CAD.smoking = CAD.group[CAD.group$smoking==1,]
plot(factor(CAD.smoking$HTN, labels = c('False', 'True')), main = 'CAD - True: HTN')
plot(factor(CAD.smoking$DM, labels = c('False', 'True')), main = 'CAD - True: DM')

par(mfrow = c(1,1))
CAD.HTN = CAD.group[CAD.group$HTN==1,]
plot(factor(CAD.HTN$DM, labels = c('False', 'True')), main = 'CAD - True: DM')

### 특이값 탐색
### 연속형 변수에 대해 RD: MCD 사용

### PCA 로 차원을 축소
pca.out = princomp(con_data[,2:9], cor = TRUE)
summary(pca.out)
### 차원 축소가 적절하게 되지 않아 생략

### Robust MCD Distance 계산
### FAST MCD 알고리즘 활용 mu & S:cov 추정
X = con_data
mcd.out = cov.mcd(X, method = 'mcd', nsamp = 'best')
mu = mcd.out$center; cov = mcd.out$cov

### RD(X) 계산 
### RD(X) = sqrt(diag((X.cent) * solve(cov) * t(X.cent)))
X.mat = matrix(stack(X)[,1], nrow = nrow(X))
X.cent = t(apply(X.mat, 1, FUN = function(x) x-mu))
rd = sqrt(diag(X.cent %*% solve(cov) %*% t(X.cent)))

### RD(X) 은 자유도 (8 - 1) 인 카이제곱 분포를 따른다.
### 이때 df =  7, p = 0.975 를 기준으로 활용
lim = sqrt(qchisq(p = 0.975, df = 7))

### RD(X) > lim 을 높은 지렛점으로 판별
HL = as.integer(rd > lim)
idx.HL = which(HL == 1)

### Mahalanobis Distance (MD) 계산 
mu.MD = matrix(apply(X.mat, 2, mean)); cov.MD = cov(X.mat)
X.cent.MD = t(apply(X.mat, 1, FUN = function(x) x-mu.MD))
md = sqrt(diag(X.cent.MD %*% solve(cov.MD) %*% t(X.cent.MD)))

HL.MD = as.integer(md > lim)
idx.HL.MD = which(HL.MD == 1)

### RD 와 MD 비교 시각화
plot(rd ~ md, ylab = 'MCD', xlab = 'MD', main = 'MCD vs RD',
     cex = 0.8, col = 'skyblue')
abline(h = lim, v = lim, col = 'red')
identify(rd ~ md, cex = 0.8)


### Robust 한 높은 지렛점 시각화
### 시각화를 위해 PCA 를 활용해 2차원으로 축소
### PCA 는 단위 문제 때문에 상관행렬을 활용
X.pca = data.frame(princomp(X)$scores[,1:2])
names(X.pca) = c('V1', 'V2')
X.pca['HL'] = HL
plot(V1 ~ V2, cex = 0.5, col = HL+1, pch = HL, data = X.pca)
### 차원 축소가 제대로 되지 않아 특징적으로 나뉘지 않는다.




