library(dplyr)
library(lubridate)

# Load dataset
df <- read.csv("df.csv", na = "null")

pen_surface <- read.csv("pen_surface.csv") %>%
    mutate(
        DateMin = as.Date(DateMin, format = "%Y-%m-%d"),
        DateMax = as.Date(DateMax, format = "%Y-%m-%d"),
    ) %>%
    select(-X)


#Already filtered: Patio == 0, Endflocksize < 80000, Mortality <= 15%, >=2017,
#NumberOfFlockIDs == 1, TimesThinned <3, 1 slaughter day, no missing slaughter weight
# slow-growing flocks are still included

# First dataset re-ordening (for all following datasets)
df_cleaned <- df %>%
    mutate(
        #remove 0 scores (but don't filter flocks)
        FootpadLesions = na_if(FootpadLesions, 0),
        FootpadLesionsLargestBatch = na_if(FootpadLesionsLargestBatch, 0),
        Farm = factor(FarmIdentification),
        Farmhouse = paste(FarmIdentification, Pen, sep = "_"),
        SlaughterDate = as.Date(SlaughterDate, format = "%Y-%m-%d"),
        SlaughterMonth = as.numeric(format(SlaughterDate, format= "%m")),
        SlaughterYear = as.numeric(format(SlaughterDate, format = "%Y")),
        StartFlockSize = as.integer(round(EndFlockSize/(1-(Mortality/100)), 0)),
        TotalDead = as.integer(StartFlockSize - EndFlockSize),
        
        HatchMonth = as.numeric(substring(HatchDate, 6, 7)),
        HatchDate = as_date(ymd_hms(HatchDate)),
        SlaughterYear = as.numeric(substring(SlaughterDate, 6, 7)),
        
        Farmhouse = factor(paste(Farm, Pen, sep = "_")),
        Farm = factor(Farm),
        Pen = factor(Pen),
        VetId = factor(VetId),
        HatchYear = factor(HatchYear),
        HatchQuarter = factor(HatchQuarter),
        HatchMonth = factor(HatchMonth),
        
        PlacementDiff = NumberPlaced - StartFlockSize
    ) %>%
    rename(
        TotalNProcessed = EndFlockSize,
        EndFlockSize = NumberSlaughtered
    ) %>%
    filter(Type != "free range and organic") %>%
    select(-FarmIdentification) %>%
    group_by(Farm, Pen) %>%
    arrange(Flock) %>%
    mutate(
        PrevTimesThinned = lag(TimesThinned),
        PrevThinned = ifelse(PrevTimesThinned == 0, 0, 1)
    ) %>%
    ungroup()

# calculate yearly FPL scores before filtering out slow-growing
dfPrevYear <- df_cleaned %>%
    group_by(Farm, Pen, SlaughterYear) %>%
    summarise(
        YearlyFPL = mean(FootpadLesionsLargestBatch, , na.rm = TRUE),
        YearlyFPLGroup = ifelse(
            YearlyFPL <= 80, 
            'Good', 
            ifelse(
                YearlyFPL <= 120, 
                'High', 
                'VeryHigh')
            ),
        .groups = 'drop'
    ) %>%
    group_by(Farm, Pen) %>%
    arrange(SlaughterYear, .by_group = TRUE) %>%
    mutate(
        PrevYearlyFPL = lag(YearlyFPL),
        PrevYearlyFPLGroup = ifelse(
            PrevYearlyFPL <= 80, 
            'Good', 
            ifelse(
                PrevYearlyFPL <= 120, 
                'High', 
                'VeryHigh')
            )
    )

# check number of flocks before additional filters on conv flocks
# and calculate quartiles
df_conv_unfiltered <- df_cleaned %>%
    filter(
        Type == "conventional"
    ) %>%
    left_join(dfPrevYear) %>%
    #add surface from other dataset
    left_join(pen_surface,
              by = join_by(HatchDate >= DateMin, HatchDate < DateMax,
                           KIPNumber, Pen)
    ) %>%
    #calculate stocking dens from gross surface
    mutate(
        FinalStockDens = SlaughterWeight / GrossSurface,
        PlacementStockDens = StartFlockSize / GrossSurface,
    ) %>%
    filter(
        !is.na(PlacementStockDens)
    )

nrow(df_conv_unfiltered) #31385

hist(df_conv_unfiltered$AgeAtSlaughter)
quantile(df_conv_unfiltered$AgeAtSlaughter, 0.01)

df_conv_2 <- df_conv_unfiltered %>%
    filter(
        AgeAtSlaughter > 34
    )

nrow(df_conv_2)

# Quantiles used for filter (based on df_conv_unfiltered)
quantile(df_conv_2$FinalStockDens, 0.01)
quantile(df_conv_2$FinalStockDens, 0.99)

quantile(df_conv_2$PlacementStockDens, 0.01)
quantile(df_conv_2$PlacementStockDens, 0.99)


# Only conventional, including not-thinned flocks.
AnalysisData_conv <- df_conv_2 %>%
    filter(
        FinalStockDens > quantile(df_conv_2$FinalStockDens, 0.01),
        FinalStockDens < quantile(df_conv_2$FinalStockDens, 0.99),
        PlacementStockDens > quantile(df_conv_2$PlacementStockDens, 0.01),
        PlacementStockDens < quantile(df_conv_2$PlacementStockDens, 0.99),
        
        !is.na(VetId)
    )

nrow(AnalysisData_conv)
summary(AnalysisData_conv$FinalStockDens)


# Save df

df <- AnalysisData_conv %>%
    select(
        Farm,
        KIPNumber,
        Pen,
        Flock,
        VetId,
        NumberOfHouses,
        GrossSurface,
        BuildYear,
        PlacementStockDens,
        FinalStockDens,
        HatchYear,
        HatchQuarter,
        TimesThinned,
        FractionThinned,
        AgeAtSlaughter,
        FirstThinningAge,
        SecondThinningAge,
        SlaughterWeightPerChicken,
        WeightAtFirstThinningPerChicken,
        AntibioticsWeek1,
        AntibioticsAfterWeek1,
        PrevYearlyFPLGroup,
        YearlyFPLGroup,
        MortalityAtFirstThinning,
        Mortality,
        StartFlockSize,
        FootpadLesions,
        PrevTimesThinned,
        PrevThinned,
        NumberPlaced,
        PlacementDiff
    )

write.csv(df, "data_conv.csv", row.names = FALSE)

# How many  missing values in outcomes?
nrow(df)
colSums(is.na(df))


