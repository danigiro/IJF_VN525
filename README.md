
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IJF_VN525

<!-- badges: start -->
<!-- badges: end -->
<!-- badges: start -->
<!-- badges: end -->

The monthly time series reconciliation of the Australian Tourism Flows
disaggregated by geographic divisions and purpose of travel, firstly
studied by Wickramasuriya et al. (2019), and used by Hollyman et
al. (2021)

## Files

-   `VN525.RData`: Australian Tourism Flows dataset
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
        -   `VN555_cslccd_bLCC.R`: LxCC with ETS bts (x = number of
            level)
        -   `VN555_cslccd_mean_red.R`: LCC with SA bts
        -   `VN555_cslccd_mean.R`: CCC with SA bts
        -   `VN555_cslccd_mLCC.R`: LxCC with SA bts (x = number of
            level)
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

Hollyman, R., Petropoulos, F., Tipping, M.E., 2021. Understanding
Forecast Reconciliation. European Journal of Operational Research 294,
149–160. <doi:10.1016/j.ejor.2021.01.017>.

Wickramasuriya, S.L., Athanasopoulos, G., Hyndman, R.J., 2019. Optimal
Forecast Reconciliation for Hierarchical and Grouped Time Series Through
Trace Minimization. Journal of the American Statistical Association 114,
804–819. <doi:10.1080/01621459.2018.1448825>.
