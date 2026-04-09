
## NMI bar plot
```R
######_----NMI

library(tidyverse)
df <- data.frame(
  method = rep(c("HVG","RBP","TF","RBP_TF"), each=10),
  NMI = c(
    # HVG
    0.8675, 0.8672, 0.8665, 0.8678, 0.8658,
    0.8669, 0.8684, 0.8660, 0.8675, 0.8671,
    
    # RBP
    0.9767, 0.9222, 0.9757, 0.9761, 0.9165,
    0.9628, 0.9805, 0.9490, 0.9456, 0.9746,
    
    # TF
    0.9861, 0.9536, 0.9759, 0.9135, 0.9772,
    0.9870, 0.9150, 0.9188, 0.9830, 0.9834,
    
    # RBP + TF
    0.9894, 0.9816, 0.9830, 0.9883, 0.9865,
    0.9889, 0.9877, 0.9880, 0.9885, 0.9892
  )
)
# method 顺序
df$method <- factor(df$method,
                    levels=c("HVG","RBP","TF","RBP_TF"))

# =========================
# 2. 计算 mean + sd
# =========================
summary_df <- df %>%
  group_by(method) %>%
  summarise(
    mean = mean(NMI),
    sd   = sd(NMI)
  )

# =========================
# 3. 作图
# =========================
ggplot(summary_df, aes(x=method, y=mean, fill=method)) +
  
  geom_bar(
    stat="identity",
    width=0.6,
    color="black"
  ) +
  
  geom_errorbar(
    aes(ymin=mean-sd, ymax=mean+sd),
    width=0.15,
    size=0.8
  ) +
  
  geom_jitter(
    data=df,
    aes(x=method, y=NMI),
    width=0.05,
    size=2.5,
    color="black"
  ) +
  
  scale_fill_manual(values=c(
    "HVG"="#4DBBD5",
    "TF"="#E64B35",
    "RBP"="#00A087",
    "RBP_TF"="#3C5488"
  )) +
  
  theme_classic() +
  
  labs(
    x="Feature representation",
    y="Normalized Mutual Information (NMI)"
  ) +
  
  theme(
    legend.position="none",
    axis.text=element_text(size=12),
    axis.title=element_text(size=13)
  )

```

### Wilcox test
```R
#----# Data (NMI)
df_nmi <- data.frame(
  replicate = c("seed1231","seed1232","seed1233","seed1234","seed1235",
                "seed1236","seed1237","seed1238","seed1239","seed12310"),
  
  HVG = c(0.8675, 0.8672, 0.8665, 0.8678, 0.8658,
          0.8669, 0.8684, 0.8660, 0.8675, 0.8671),
  
  RBP = c(0.9767, 0.9222, 0.9757, 0.9761, 0.9165,
          0.9628, 0.9805, 0.9490, 0.9456, 0.9746),
  
  TF = c(0.9861, 0.9536, 0.9759, 0.9135, 0.9772,
         0.9870, 0.9150, 0.9188, 0.9830, 0.9834),
  
  RBP_TF = c(0.9894, 0.9816, 0.9830, 0.9883, 0.9865,
             0.9889, 0.9877, 0.9880, 0.9885, 0.9892)
)

#Friedman test
mat_nmi <- as.matrix(df_nmi[, -1])
friedman.test(mat_nmi)


#paired Wilcoxon（two-sided）
wilcox.test(df_nmi$RBP, df_nmi$HVG, paired = TRUE, exact = TRUE)
wilcox.test(df_nmi$TF, df_nmi$HVG, paired = TRUE, exact = TRUE)
wilcox.test(df_nmi$RBP_TF, df_nmi$HVG, paired = TRUE, exact = TRUE)

wilcox.test(df_nmi$RBP, df_nmi$TF, paired = TRUE, exact = TRUE)
wilcox.test(df_nmi$RBP_TF, df_nmi$RBP, paired = TRUE, exact = TRUE)
wilcox.test(df_nmi$RBP_TF, df_nmi$TF, paired = TRUE, exact = TRUE)

````




