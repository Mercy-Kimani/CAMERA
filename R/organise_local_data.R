#' CAMERA_local class
#'
#' @description
#' A simple wrapper function for importing data from local files for use with the CAMERA class.
#' @export
CAMERA_local <- R6::R6Class("CAMERA_local", list(
    metadata = NULL,
    ld_ref = NULL,
    mc.cores = NULL,
    plink_bin = NULL,
    radius = NULL,
    minmaf = NULL,
    pthresh = NULL,
    instrument_raw = NULL,
    instrument_outcome = NULL,
    instrument_regions = NULL,
    instrument_outcome_regions = NULL,

    # Methods
    initialize = function(metadata, ld_ref, plink_bin, mc.cores=1, radius = 25000, pthresh = 5e-8, minmaf=0.01) {
        self$metadata <- metadata
        self$plink_bin <- plink_bin
        self$radius <- radius
        self$mc.cores <- mc.cores
        self$pthresh <- pthresh
        self$minmaf <- minmaf
        self$ld_ref <- ld_ref
    },

    standardise = function(d, ea_col="ea", oa_col="oa", beta_col="beta", eaf_col="eaf", chr_col="chr", pos_col="pos", vid_col="vid") {
        toflip <- d[[ea_col]] > d[[oa_col]]
        d[[eaf_col]][toflip] <- 1 - d[[eaf_col]][toflip]
        d[[beta_col]][toflip] <- d[[beta_col]][toflip] * -1
        temp <- d[[oa_col]][toflip]
        d[[oa_col]][toflip] <- d[[ea_col]][toflip]
        d[[ea_col]][toflip] <- temp
        d[[vid_col]] <- paste0(d[[chr_col]], ":", d[[pos_col]], "_", d[[ea_col]], "_", d[[oa_col]])
        d
    },

    read_file = function(m, minmaf=0.01) {
        stopifnot(nrow(m) == 1)
        stopifnot(file.exists(m$fn))
        a <- data.table::fread(m$fn)
        message("Read ", nrow(a), " rows")
        b <- tibble(
            chr = a[[m$chr_col]],
            pos = as.numeric(a[[m$pos_col]]),
            eaf = as.numeric(a[[m$eaf_col]]),
            beta = as.numeric(a[[m$beta_col]]),
            se = as.numeric(a[[m$se_col]]),
            pval = as.numeric(a[[m$pval_col]]),
            ea = a[[m$ea_col]],
            oa = a[[m$oa_col]]
        ) %>% 
        dplyr::filter(eaf > minmaf & eaf < (1-minmaf)) %>%
        self$standardise()
        return(b)
    },

    pool_tophits = function(rawdat, tophits, metadata, radius = 250000, pthresh = 5e-8, mc.cores = 10) {
        regions <- GRanges(
            seqnames = tophits$chr,
            ranges = IRanges(start=tophits$pos-radius, end=tophits$pos+radius),
            vid=tophits$vid, 
            pop=tophits$pop,
            trait=tophits$trait
        )
        region_list <- lapply(unique(tophits$trait), function(tr) {
            temp <- reduce(subset(regions, trait == tr))
            temp$trait <- tr
            temp
        })

        # Extract regions
        region_extract <- lapply(1:length(region_list), function(tr) {
            region <- region_list[[tr]]
            lapply(1:length(region), function(i) {
                message(i, " of ", length(region))
                lapply(1:nrow(metadata), function(j) {
                    subset(rawdat[[j]], chr == as.character(seqnames(region)[i]) & pos <= end(region)[i] & pos >= start(region)[i]) %>%
                        mutate(trait = metadata$trait[j], pop = metadata$pop[j], id = metadata$id[j])
                }) %>% dplyr::bind_rows()
            })
        })

        pool <- lapply(1:length(region_extract), function(tr) {
            region <- region_extract[[tr]]
            lapply(1:length(region), function(i) {
                target_trait <- region_list[[tr]]$trait[1]
                a <- region[[i]]
                k <- a %>% dplyr::group_by(vid) %>% dplyr::summarise(nstudies=n(), .groups="drop")
                a <- dplyr::left_join(a, k, by="vid")
                k <- a %>% dplyr::filter(trait == target_trait) %>%
                    dplyr::group_by(nstudies) %>% 
                    dplyr::summarise(minp = min(pval), .groups="drop") %>% 
                    dplyr::filter(minp < pthresh)
                a <- subset(a, nstudies %in% k$nstudies)
                k <- subset(a, trait == target_trait) %>% 
                    mutate(z = abs(beta)/se) %>%
                    {subset(., z==max(z))$vid[1]}
                a <- subset(a, vid == k) %>% mutate(target_trait=target_trait)
                return(a)
            }) %>% dplyr::bind_rows()
        }) %>% dplyr::bind_rows()

        region_list <- lapply(region_list, as_tibble)
        list(region_list=region_list, region_extract=region_extract, tophit_pool=pool)
    },

    organise_data = function(metadata=self$metadata, plink_bin=self$plink_bin, ld_ref=self$ld_ref, pthresh=self$pthresh, minmaf = self$minmaf, radius = self$radius, mc.cores = self$mc.cores, rawdat = NULL) {
        if(is.null(rawdat)) {
            rawdat <- lapply(1:nrow(metadata), function(i) self$read_file(metadata[i,]))
        }

        tophits <- lapply(1:nrow(metadata), function(i) {
            x <- rawdat[[i]] %>% dplyr::filter(pval < pthresh) %>% dplyr::mutate(rsid = vid)
            if(nrow(x) > 1) {
                ieugwasr::ld_clump(x, plink_bin=plink_bin, bfile=subset(ld_ref, pop == metadata$pop[i])$bfile) %>%
                    dplyr::select(-rsid) %>%
                    dplyr::mutate(pop=metadata$pop[i], trait=metadata$trait[i])
            } else {
                NULL
            }
        }) %>% dplyr::bind_rows()

        self$pool_tophits(rawdat, tophits, metadata, radius = radius, pthresh = pthresh, mc.cores = mc.cores)
    },

    fixed_effects_meta_analysis_fast = function(beta_mat, se_mat) {
        w <- 1 / se_mat^2
        beta <- rowSums(beta_mat * w, na.rm=TRUE) / rowSums(w, na.rm=TRUE)
        se <- sqrt(1 / rowSums(w, na.rm=TRUE))
        pval <- pnorm(abs(beta / se), lower.tail = FALSE)
        return(pval)
    },

    organise = function() {
        metadata <- self$metadata
        ld_ref <- self$ld_ref

        metadata_exp <- metadata[metadata$what == "exposure", ]
        metadata_out <- metadata[metadata$what == "outcome", ]

        exposure_trait <- metadata_exp$trait[1]
        outcome_trait <- metadata_out$trait[1]

        # Read exposure datasets
        rawdat <- lapply(1:nrow(metadata_exp), function(i) self$read_file(metadata_exp[i, ]))

        # Process each outcome trait
        outcome_traits <- unique(metadata_out$trait)
        out <- lapply(outcome_traits, function(trt) {
            message("Processing outcome: ", trt)
            temp <- metadata_out[metadata_out$trait == trt, ]
            rawdat_this <- lapply(1:nrow(temp), function(i) self$read_file(temp[i, ]))
            rawdat_all <- c(rawdat, rawdat_this)

            self$organise_data(subset(metadata, trait %in% c(exposure_trait, trt)),
                               self$plink_bin,
                               self$ld_ref,
                               self$pthresh,
                               self$minmaf,
                               self$radius,
                               self$mc.cores,
                               rawdat = rawdat_all)
        })

        o <- out[[1]]

        inst <- unique(subset(o$tophit_pool, target_trait == exposure_trait)$vid)
        inst_o <- unique(subset(o$tophit_pool, target_trait == outcome_trait)$vid)

        names(o$region_extract[[1]]) <- inst
        names(o$region_extract[[2]]) <- inst_o

        instrument_raw <- o$tophit_pool %>%
            dplyr::filter(target_trait == exposure_trait & trait == exposure_trait) %>%
            dplyr::rename(position = pos, nea = oa, p = pval, rsid = vid)

        instrument_outcome <- o$tophit_pool %>%
            dplyr::filter(trait == outcome_trait & target_trait == exposure_trait & vid %in% instrument_raw$rsid) %>%
            dplyr::rename(position = pos, nea = oa, p = pval, rsid = vid)

        instrument_regions <- lapply(unique(instrument_raw$rsid), function(x) {
            a <- o$region_extract[[1]][[x]] %>%
                dplyr::filter(trait == exposure_trait) %>%
                dplyr::rename(position = pos, nea = oa, p = pval, rsid = vid) %>%
                dplyr::group_by(pop) %>% dplyr::group_split() %>% as.list()
            names(a) <- sapply(a, function(z) z$id[1])
            a
        })

        instrument_outcome_regions <- lapply(unique(instrument_raw$rsid), function(x) {
            a <- o$region_extract[[1]][[x]] %>%
                dplyr::filter(trait == outcome_trait) %>%
                dplyr::rename(position = pos, nea = oa, p = pval, rsid = vid) %>%
                dplyr::group_by(pop) %>% dplyr::group_split() %>% as.list()
            names(a) <- sapply(a, function(z) z$id[1])
            a
        })

        names(instrument_regions) <- unique(instrument_raw$rsid)
        names(instrument_outcome_regions) <- unique(instrument_raw$rsid)

        self$instrument_regions <- instrument_regions
        self$instrument_outcome_regions <- instrument_outcome_regions
        self$instrument_raw <- instrument_raw
        self$instrument_outcome <- instrument_outcome
    }
))
