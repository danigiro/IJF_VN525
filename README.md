
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IJF_VN525

<!-- badges: start -->
<!-- badges: end -->

Forecast combination based reconciliation of the monthly Australian
Tourism Demand series disaggregated by geographic divisions and purpose
of travel (Di Fonzo and Girolimetto, 2022)

The forecast reconciliation is performed using **FoReco 0.2.2**.

## Files

-   `VN525.RData`: Australian Tourism Flows dataset (Visitor Nights)
-   `Ctools555.RData`: tools to balance the C matrix (VN525 -> VN555)
-   **BaseForecasts**:
    -   `ETSlev/VN525_base_ets_lev.R`: ETS base forecasts
    -   `ETSlev/VN555/VN555_base.R`: VN555 ETS levels and SA base
        forecasts
    -   `VN525_features_residuals.R`: Seasonal Average residuals
    -   `VN525_features.R`: Seasonal Average forecasts
-   **Reconciliation**:
    -   `ETSlev`:
        -   `VN555_cslccd_bCCCred.R`: LCC with ETS bts
        -   `VN555_cslccd_bLCC.R`: LxCC with ETS bts (x = level number)
        -   `VN555_cslccd_mean_red.R`: LCC with SA bts
        -   `VN555_cslccd_mean.R`: CCC with SA bts
        -   `VN555_cslccd_mLCC.R`: LxCC with SA bts (x = level number)
        -   `VN555_cslccd.R`: CCC with ETS bts and CCCH
        -   `VN525_htsrec_mean.R`: cross-sectional reconciliation with
            SA bts
        -   `VN525_htsrec.R`: cross-sectional reconciliation with ETS
            bts
        -   `VN555_CCCmix_red.R`: sample average of LCC (SA + ETS bts)
        -   `VN555_CCCmix.R`: sample average of CCC (SA + ETS bts)
        -   `VN525_csmix.R`: sample average of the cross-sectional
            reconciliation (SA + ETS bts)
    -   `function`: functions for `score/VN525_scores_ETSlev.R`
    -   `score/VN525_scores_ETSlev.R`: returns the scores’ dataset

## References

Di Fonzo, T. and Girolimetto D. (2022). *Forecast combination based
forecast reconciliation: insights and extensions*, preprint submitted to
the International Journal of Forecasting
