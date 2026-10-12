test_that("svcLMexact with a varying intercept only matches spLMexact", {
  dat <- simSpatial[1:60, ]
  co <- as.matrix(dat[, c("s1", "s2")])
  for(cf in c("exponential", "matern")){
    spp <- if(cf == "matern") list(phi = 4, nu = 1.5) else list(phi = 4)
    m1 <- spLMexact(y_gauss ~ x1 + x2, data = dat, coords = co, cor.fn = cf,
                    spParams = spp, noise_sp_ratio = 0.3, n.samples = 10,
                    loopd = TRUE, loopd.method = "exact", verbose = FALSE)
    m2 <- svcLMexact(y_gauss ~ x1 + x2 + (1), data = dat, coords = co, cor.fn = cf,
                     process.type = "independent.shared", spParams = spp,
                     noise_sp_ratio = 0.3, n.samples = 10, loopd = TRUE,
                     verbose = FALSE)
    expect_equal(m2$loopd, m1$loopd, tolerance = 1e-10)
  }
})

test_that("svcLMexact (exponential) matches stvcLMexact at a constant time", {
  dat <- simSpatial[1:60, ]
  co <- as.matrix(dat[, c("s1", "s2")])
  set.seed(1)
  m1 <- svcLMexact(y_svc ~ x1 + x2 + (x1), data = dat, coords = co,
                   cor.fn = "exponential", process.type = "independent",
                   spParams = list(phi = c(3, 6)), noise_sp_ratio = c(0.3, 1),
                   n.samples = 20, loopd = TRUE, verbose = FALSE)
  set.seed(1)
  m2 <- stvcLMexact(y_svc ~ x1 + x2 + (x1), data = dat, sp_coords = co,
                    time_coords = matrix(0, nrow(dat), 1), cor.fn = "gneiting-decay",
                    process.type = "independent",
                    sptParams = list(phi_s = c(3, 6), phi_t = c(1, 1)),
                    noise_sp_ratio = c(0.3, 1), n.samples = 20, loopd = TRUE,
                    verbose = FALSE)
  expect_equal(m1$loopd, m2$loopd, tolerance = 1e-10)
  expect_equal(m1$samples$z, m2$samples$z, tolerance = 1e-10)
})
