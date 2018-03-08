#' Calibration Analysis of Survival Predictions
#'
#' @rdname survCalibCurve
#'
#' @param data A dataset
#' @param time unquoted name of the time variable
#' @param event unquoted name of the event variable (0/1)
#' @param p unquoted name of the variable containing predicted probabilities (at
#'   t=timepoint)
#' @param timepoint Value of 'time' where probabilities are extracted (from the
#'   Kaplan-Meier curves).
#' @param showCStat Logical. If TRUE, a concordance statistic (from
#'   \code{rms::somers2()}) with confidence limits (2 std. errors) will be displayed
#' @param p.out grid of predicted probability values at which calibration curve
#'   is estimated. If NULL, a sequence of 10 values evenly spaced along the range
#'   of inputted predicted probabilities is generated.
#' @param na.rm Logical. Remove rows with missing values?
#' @param knots knots parameter passed to \code{splines::bs()}
#' @param bs.degree degree parameter passed to \code{splines::bs()}
#'
#' @details The predicted probabilities are transform via the logit prior to
#'   spline curve fitting, in order to stabilize the resulting calibration
#'   curve estimates.  Therefore, if specifying the \code{knots=} parameter to
#'   \code{bs.params}, the values must be expressed on the logit scale. (It's
#'   usually easier/sufficient to just supply \code{df} instead of \code{knots}.)
#'
#'   Spline fits in Cox proportional hazards model are prone to minor convergence
#'   issues, which impact the validity of Wald-based standard errors, tests, and
#'   confidence intervals. Hence this may produce a warning of the form "Loglik
#'   converged before variable X; beta may be infinite". See
#'   https://stat.ethz.ch/pipermail/r-help/2008-September/174201.html for details.
#'
#'   The user is advised to choose the knots for the splines very carefully,
#'   especially in the presence of highly skewed distributions of model-predicted
#'   probabilities. The default values are just a suggestion for a typical
#'   scenario where predictions are highly skewed to the right.
#'
#' @return
#' @export
#'
#' @examples

survCalibCurve <- function(data, time, event, p, timepoint,
                           showCStat=TRUE, p.out = NULL, na.rm=FALSE,
                           bs.knots=c(.01,.025,.05,.1,.2,.3), bs.degree=3) {

  time  <- enquo(time)
  event <- enquo(event)
  p     <- enquo(p)
  gv    <- group_vars(data)

  md <- data %>%
    select(group_vars(data), !!time, !!event, !!p) %>%
    rename(time=!!time, event=!!event, p=!!p) %>%
    mutate(risk.gp = as.character(cut(p, c(0,bs.knots,1), include.lowest = TRUE)))

  if(na.rm){
    md <- md %>% filter(!is.na(time) & !is.na(event) & !is.na(p))
  }
  if(any(apply(md,2,function(x) sum(is.na(x)))>0)) stop("Input data contains NAs")
  if(min(md$p)<0 | max(md$p) > 1) stop("Parameter 'p' must contain probabilities.")

  if(is.null(p.out)) {
    p.out <- seq(min(md$p), max(md$p), length.out=100)
  } else{
    if(min(p.out)<0 | max(p.out) > 1) {
      stop("Parameter 'p.out' must contain probabilities.")
    }
  }

  # this subroutine fits a spline model to get the calibration curve
  getCalibCurve <- function(data){
    g   <- survival::coxph(survival::Surv(time, event) ~
                             splines::bs(p, knots = bs.knots, degree = bs.degree),
                           data = data)
    pg  <- survival::survfit(g, newdata = tibble(p = p.out))
    tix <- which(pg$time >= timepoint)[1]
    tibble(p         = p.out,
           incidence = 1-pg$surv[tix,],
           lcl       = 1-pg$upper[tix,],
           ucl       = 1-pg$lower[tix,])
  }

  # this subroutine splits the data into risk groups and does the same
  getRiskGpCalibData <- function(data){
    rgps <- unique(data$risk.gp)
    mean.p <- group_by(data, risk.gp) %>%
      summarize(mean.p = mean(p)) %>%
      select(risk.gp, mean.p)

    h   <- survival::coxph(survival::Surv(time, event) ~ risk.gp, data = data)
    ph  <- survival::survfit(h, newdata = tibble(risk.gp = rgps))
    tix <- which(ph$time >= timepoint)[1]

    left_join(mean.p,
              tibble::tibble(risk.gp   = rgps,
                             incidence = 1-ph$surv[tix,],
                             lcl       = 1-ph$upper[tix,],
                             ucl       = 1-ph$lower[tix,]),
              by = "risk.gp")
  }

  list(
    cal = plyr::ddply(md, gv, getCalibCurve),
    risk.gp.cal = plyr::ddply(md, gv, getRiskGpCalibData),
    knots = bs.knots
  )
}


#d$pcerm <- d$pcermpct / 100; d$pcermpct <- NULL
# z <- group_by(d, sex, diabetes) %>%
#   survCalibCurve(event_weeks, ascvd_event, pcerm, 260, na.rm=TRUE)
