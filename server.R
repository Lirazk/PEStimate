# Maybe can move the withProgress to outside of the first {
library(shiny)
# Maybe I should switch to tinyplot? should be a bit more lightweight
library(ggplot2)
library(ggdark)

source("EmbryoSelection.R")

print_family_history <- function(temp) {
  cat(
    sprintf(
      "<p><u>Risk <strong>without</strong> embryo selection</u>: <strong>%.2f%% (%.2f)</strong>\n</p>",
      100*temp[1], 100*temp[6]
    )
  )
  cat(
    sprintf(
      "<p><u>Risk <strong>with</strong> embryo selection</u>: <strong>%.2f%% (%.2f)</strong>\n</p>",
      100*temp[2],
      100*temp[5]
    )
  )
  
  # cat(sprintf("<p><u>Relative risk reduction</u>: <strong>%.2f (%.2f) </strong>\n</p>", temp[3], temp[5] / temp[1]))
  cat(sprintf("<p><u>Relative risk reduction</u>: <strong>%.2f%% (%.2f) </strong>\n</p>", 100*temp[3], 100*sqrt(temp[5]^2 * temp[2]^2 / temp[1]^4 + temp[6]^2 / temp[1]^2)))
  cat(sprintf("<p><u>Absolute risk reduction</u>: <strong>%.2f%% (%.2f)</strong>\n</p>", 100*temp[4], 100*sqrt(temp[5]^2+temp[6]^2)))
  
  # cat(sprintf("<p><u>Couples needed to screen</u>: <strong>%.0f (%.4f)</strong>\n</p>", ceiling(1/temp[4]), temp[5] / temp[4]^2))
  cat(sprintf("<p><u>Couples needed to screen</u>: <strong>%.0f (%.4f)</strong>\n</p>", ceiling(1/temp[4]), sqrt((temp[5]^2 + temp[6]^2) / temp[4]^4)))
  
  # cat(sprintf("<p style = \"color:red\">Based on %d draws. Estimated standard deviation in parentheses.</p>", input$samples))
  cat(sprintf("<p style = \"color:red\">Based on simulation. Estimated standard deviation in parentheses.</p>"))
  if(temp[1] < temp[2]) {
    cat(sprintf("<b><p style = \"color:red\">Baseline risk is smaller than the strategy risk. Either the sample size is too small, or the baseline and strategy risk are almost identical.</p></b>"))
  }
}

print_result <- function(temp, K) {
  cat(
    sprintf(
      "<p><u>Risk <strong>without</strong> embryo selection</u>: <strong>%.2f%%</strong>\n</p>",
      # K/100
      K
    )
  )
  cat(
    sprintf(
      "<p><u>Risk <strong>with</strong> embryo selection</u>: <strong>%.2f%%</strong>\n</p>",
      # K/100 - temp * K/100
      K - temp * K
    )
  )
  cat(sprintf("<p><u>Relative risk reduction</u>: <strong>%.2f%%</strong>\n</p>", 100*temp))
  cat(sprintf("<p><u>Absolute risk reduction</u>: <strong>%.2f%%</strong></p>", temp * K))
  
  # TODO: is it correct?
  cat(sprintf("<p><u>Couples needed to screen</u>: <strong>%.0f</strong></p>", ceiling(100/(temp * K))))
}

print_result_conditional <- function(temp) {
  cat(
    sprintf(
      "<p><u>Risk <strong>without</strong> embryo selection</u>: <strong>%.2f%%</strong>\n</p>",
      100*temp$baseline
    )
  )
  cat(
    sprintf(
      "<p><u>Risk <strong>with</strong> embryo selection</u>: <strong>%.2f%%</strong>\n</p>",
      100*temp$risk
    )
  )
  cat(sprintf("<p><u>Relative risk reduction</u>: <strong>%.2f%%</strong>\n</p>", 100*temp$rr))
  cat(sprintf("<p><u>Absolute risk reduction</u>: <strong>%.2f%%</strong></p>", 100*temp$rr * temp$baseline))
  cat(sprintf("<p><u>Couples needed to screen</u>: <strong>%.0f</strong></p>", ceiling(1/(temp$rr * temp$baseline))))
}

server <- function(input, output, session) {
  accepeted <- F
  modal <- modalDialog(h1("Important notes:"), 
                       HTML("<ol>
       <li>Our risk estimates are based on the liability threshold model from statistical genetics. <strong>They are not directly based on real epidemiological data</strong>.</li>
       <li>We assume that patients screen their embryos for a <strong>single disease</strong>.</li>
       <li>Under the model, <strong>a disease is assumed to have an underlying, continuous liability</strong>. The liability is normally distributed and it represents the sum of additive genetic and non-genetic risk factors. Under the model, individuals are affected if their liability exceeds a threshold.</li>
       <li>The model assumes that the <strong>disease status is binary</strong>: an adult individual is either affected or unaffected. Both the prevalence parameter and the estimated risk represent the proportion of affected adults.</li>
       <li>The model assumes that a PRS explains a proportion R² of the variance in the liability in the population. This parameter, which quantifies the accuracy of the PRS, <strong>should reflect the accuracy of risk prediction between siblings for the setting of interest</strong>. Factors that may reduce accuracy compared to values reported in the literature include, for example, differences in sequencing platforms, genotyping errors, application in populations of non-European descent, population structure, parental environment correlated with child’s genotype, and lower accuracy of the PRS in future years.</li>
       <li>All embryos are assumed to be <strong>euploid blastocysts</strong>. We assume a single embryo transfer. The outcome of the transfer is a live birth with either probability 1 or a prespecified probability. The <strong>live birth probability is assumed to be identical between embryos and patients</strong>. In the models with a random number of births, it is assumed that <strong>at least one birth was achieved</strong>.</li>
       <li>The app does not provide information on <strong>the likelihood of increasing the risk of other diseases due to selection</strong>.</li>
       <li><strong>Other assumptions</strong> underlying the model are discussed in our <a href=\"https://doi.org/10.7554/elife.64716\" target=\"_blank\" rel=\"noopener noreferrer\">paper</a>.</li>
       <li>Finally, note that screening IVF embryos for polygenic risk scores is associated with <strong>multiple ethical, social, and legal concerns</strong>. For a comprehensive review of epidemiological, clinical, and ethical considerations, please see our <a href=\"https://doi.org/10.1093/humupd/dmae012\" target=\"_blank\" rel=\"noopener noreferrer\">review</a>.</li>
       </ol>"),
                       p("Please confirm that:"),
                       HTML("<ol>
       <li>You read and understood the above notes.</li>
       <li>You understand that the application is intended for research purposes only and is not intended to guide clinical decision making.</li>
       </ol>"), checkboxInput("accept", "Accept"), easyClose = F, footer = NULL, size = "l")

  observeEvent(input$accept, {
    req(input$accept)
    removeModal()
    accepeted <<- T
  })
  
  observeEvent(input$type2, {
    if(input$type2 == 'Conditional') {
      shinyjs::enable("qf2")
      shinyjs::enable("qm2")
    }
    else {
      shinyjs::disable("qf2")
      shinyjs::disable("qm2")
    }
  })
  
  observeEvent(input$lowestexclude2, {
    if(input$lowestexclude2 == 'Lowest') {
      shinyjs::disable("q2")
    }
    else {
      shinyjs::enable("q2")
    }
  })
  
  values <- reactiveValues(disable = FALSE)
  
  observeEvent(input$parents_n | input$siblings_n, {
    if (input$parents_n == 0 & input$siblings_n == 0) {
      shinyjs::disable("h2")
      values$disable <- T
    }
    else {
      shinyjs::enable("h2")
      if (input$r2 >= input$h2) {
        updateSliderInput(session, "r2", value = input$h2)
      }
      values$disable <- F
    }
  })
  
  # Is it still needed?
  choice_names <- c("PRS accuracy (R-squared)", "Disease prevalence", "Number of live births", 
                    "Percentile")
  choice_values <- c("r2", "Disease prevalence", "Number of embryos", "Percentile")
  
  observeEvent(input$det_random, {
    if (input$det_random == "Binomial") {
      choice_names[3] <- "Number of euploid embryos"
      updateSliderInput(session, "N", max = 20,
                            label = "Number of euploid embryos")
    }
    else if (input$det_random == "Poisson") {
      choice_names[3] <- "Expected number of euploid embryos"
      updateSliderInput(session, "N", max = 20, 
                            label = "Expected number of euploid embryos")
    }
    else {
      choice_names[3] <- "Number of live births"
      updateSliderInput(session, "N", max = 10,
                            label = "Number of live births")
    }
    updateRadioButtons(session, "x_var", choiceNames = choice_names[1:3], 
                       choiceValues = choice_values[1:3], inline = T, 
                       selected = input$x_var)
  })
  
  observeEvent(input$lowestexclude, {
    # input$det_random
    if (input$lowestexclude == "Lowest") {
      if (input$det_random == "Binomial") {
        choice_names[3] <- "Number of euploid embryos"
        updateSliderInput(session, "N", max = 20,
                          label = "Number of euploid embryos")
      }
      else if (input$det_random == "Poisson") {
        choice_names[3] <- "Expected number of euploid embryos"
        updateSliderInput(session, "N", max = 20, 
                          label = "Expected number of euploid embryos")
      }
      else {
        choice_names[3] <- "Number of live births"
        updateSliderInput(session, "N", max = 10,
                          label = "Number of live births")
      }
      updateRadioButtons(session, "x_var", choiceNames = choice_names[1:3], 
                         choiceValues = choice_values[1:3], inline = T, 
                         selected = input$x_var)
    }
    else {
      choice_names[3] <- "Number of live births"
      updateRadioButtons(session, "x_var", choiceNames = choice_names,
                                     choiceValues = choice_values, inline = T)
      updateSliderInput(session, "N", max = 10, label = "Number of live births")
    }
  })
  
  filter_input <- debounce(reactive({
    list(r2 = input$r,
         K2 = input$K,
         N = input$N,
         q = input$q,
         p = input$p_lb)
  }), 500)
  
  observeEvent(input$r2, {
    if (input$r2 >= input$h2 & !values$disable) {
      updateSliderInput(session, "r2", value = input$h2)
    }
  })
  
  observeEvent(input$parents_n | input$sick_parents, {
    # updateSliderInput(session, "sick_parents", max=input$parents_n,
                      # value=pmin(input$parents_n, input$sick_parents))
    if (input$sick_parents > input$parents_n) {
      # input$sick_parents <- input$parents_n
      updateSliderInput(session, "sick_parents", value = input$parents_n)
    }
  })
  
  observeEvent(input$siblings_n | input$sick_siblings, {
    # updateSliderInput(session, "sick_siblings", max=input$siblings_n,
                      # value=pmin(input$siblings_n, input$sick_siblings))
    if (input$sick_siblings > input$siblings_n) {
      # input$sick_siblings <- siblings_n
      updateSliderInput(session, "sick_siblings", value = input$siblings_n)
    }
  })
  
  mode_color <- reactive({
    if (input$dark_mode == "dark") {
      "white"
    }
    else {
      "black"
    }
  })
  
  dark_theme <- dark_mode(theme_minimal())
  
  plot_theme <- reactive({
    if (input$dark_mode == "dark") {
      # dark_theme_minimal() +
      dark_theme +
      # theme_dark() #+
        theme(text = element_text(color = "white"))
    } else {
      theme_minimal() +
        theme(text = element_text(color = "black"))
    }
  })

  output$distPlot <- renderPlot({withProgress({
    if (!accepeted) showModal(modal)
    
    
    filts <- filter_input()
    # output$distPlot <- renderPlotly({
    # selectInput("x_var", "Variable for x axis", choices = c("r2", "K", "N")),
    if (input$lowestexclude == "Lowest") {
      subtitle <- "Lowest PRS strategy"
    }
    else {
      subtitle <- "Exclude high PRS strategy"
    }
    
    binomial_model <- input$det_random == "Binomial"
    poisson_model <- input$det_random == "Poisson"
    
    # if(binomial_model) {
    #   updateSliderTextInput(session, "N", label = "Number of euploid embryos")
    #   # choice_names[3] <-  "Number of euploid embryos"
    # }
    # else if(poisson_model) {
    #   updateSliderTextInput(session, "N", label = "Expected number of euploid embryos")
    #   # choice_names[3] <-  "Expected number of euploid embryos"
    # }
    # else {
    #   updateSliderTextInput(session, "N", label = "Number of live births")
    #   # choice_names[3] <-  "Number of live births"
    # }

    # if (input$lowestexclude == "Lowest") {
    #   updateRadioButtons(session, "x_var", choiceNames = choice_names[1:3], inline = T)
    # }
    # else {
    #   updateRadioButtons(session, "x_var", choiceNames = choice_names, inline = T)
    # }
        
    selected_x <- input$x_var
    if(selected_x == "Number of embryos") {
      x <- 2:20
      if (input$lowestexclude == "Lowest") {
        # y <- sapply(x, function(x) risk_reduction_lowest(filts$r, filts$K, n = x))
        y <- sapply(x, function(x) risk_reduction_lowest(filts$r, filts$K/100, n = x))
        # if (input$relative_abs == "Absolute risk") y  <- y * input$K
        if (input$det_random == "Binomial") {
          # print(filts$p)
          y <- sapply(x, function(x) risk_reduction_lowest_bin(filts$r, filts$K/100, n = x, 
                                                               filts$p))
          # print(sapply(x, function(x) risk_reduction_lowest_bin(filts$r, filts$K, n = x, 
          #                                                       filts$p)))
          # print(sapply(x, function(x) risk_reduction_lowest(filts$r, filts$K/100, n = x)))
        }
        else if (input$det_random == "Poisson") {
          y <- sapply(x, function(x) risk_reduction_lowest_pois(filts$r, filts$K/100, 
                                                                x * filts$p))
        }
      }
      else {
        # y <- sapply(x, function(x) risk_reduction_exclude(r2 = filts$r, K = filts$K, q = filts$q, n = x))
        y <- sapply(x, function(x) risk_reduction_exclude(r2 = filts$r, K = filts$K/100, q = filts$q, n = x))
        # if (input$relative_abs == "Absolute risk") y  <- y * input$K
      }
            
      # subtitle <- "Lowest strategy"
      x_lab <- "Number of live births"
      if (input$lowestexclude == "Lowest") {
        if (binomial_model) {
          x_lab <- "Number of euploid embryos"
        }
        else if (poisson_model) {
          x_lab <- "Expected number of euploid embryos"
        }
      }
    }
    else if(selected_x == "Disease prevalence") {
      x <- exp(seq(log(0.001), log(0.3), length = 50))
      # y <- sapply(x, function(x) risk_reduction_lowest(input$r, x, n = input$N))
      if (input$lowestexclude == "Lowest") {
        y <- sapply(x, function(x) risk_reduction_lowest(filts$r, x, n = filts$N))
        if (input$det_random == "Binomial") {
          y <- sapply(x, function(x) risk_reduction_lowest_bin(filts$r, x, n = filts$N, input$p_lb))
        }
        else if (input$det_random == "Poisson") {
          y <- sapply(x, function(x) risk_reduction_lowest_pois(filts$r, x, filts$N*input$p_lb))
        }
        # if (input$relative_abs == "Absolute risk") y  <- y * x
      }
      else {
        y <- sapply(x, function(x) risk_reduction_exclude(r2 = filts$r, K = x, q = filts$q, n = filts$N))
        # if (input$relative_abs == "Absolute risk") y  <- y * x
      }
      # y <- binomial_random(y, input$p_lb, filts$N)
      # subtitle <- "Lowest strategy"
      x_lab <- "Disease prevalence"
    }
    else if(selected_x == "r2" || selected_x == "$$R^2$$") {
      x <- seq(0.01, 1, length = 50)
      # y <- sapply(x, function(x) risk_reduction_lowest(x, input$K, n = input$N))
      if (input$lowestexclude == "Lowest") {
        # y <- sapply(x, function(x) risk_reduction_lowest(x, filts$K, n = input$N))
        y <- sapply(x, function(x) risk_reduction_lowest(x, filts$K/100, n = input$N))
        if (input$det_random == "Binomial") {
          y <- sapply(x, function(x) risk_reduction_lowest_bin(x, filts$K/100, n = filts$N, input$p_lb))
        }
        else if (input$det_random == "Poisson") {
          y <- sapply(x, function(x) risk_reduction_lowest_pois(x, filts$K/100, filts$N*input$p_lb))
        }
        # if (input$relative_abs == "Absolute risk") y  <- y * input$K
      }
      else {
        # y <- sapply(x, function(x) risk_reduction_exclude(x, filts$K, filts$q, n = filts$N))
        y <- sapply(x, function(x) risk_reduction_exclude(x, filts$K/100, filts$q, n = filts$N))
        # if (input$relative_abs == "Absolute risk") y  <- y * input$K
      }
      # y <- binomial_random(y, input$p_lb, filts$N)
      # subtitle <- "Lowest strategy"
      # x_lab <- "PRS r^2"
      # x_lab <- expression("PRS accuracy (R^2)")
      x_lab <- expression(paste("PRS accuracy (", R^2, ")"))
    }
    else if (selected_x == "Percentile") {
      # q
      x <- seq(0.01, 1, length = 50)
      # y <- sapply(x, function(x) risk_reduction_exclude(filts$r, filts$K, x, n = filts$N))
      y <- sapply(x, function(x) risk_reduction_exclude(filts$r, filts$K/100, x, n = filts$N))
      # if (input$relative_abs == "Absolute risk") y  <- y * input$K
      # subtitle <- "Exclude strategy"
      x_lab <- "Percentile to exclude"
    }
    else {
      print("Error!")
      print(selected_x)
    }
    
    if(selected_x == "Disease prevalence") {
      ggplot() +
        geom_point(mapping = aes(x, 100*y, color = "Relative risk"), size = 3) +
        geom_line(mapping = aes(x, 100*y, color = "Relative risk")) +
        geom_point(mapping = aes(x, 100*y*x, color = "Absolute risk"), size = 3) +
        geom_line(mapping = aes(x, 100*y*x, color = "Absolute risk")) +
        # theme_minimal() +
        plot_theme() +
        labs(
          title = "Risk reduction",
          subtitle = subtitle,
          x = x_lab,
          y = "Risk reduction (%)"
        ) +
        theme(
          plot.subtitle = element_text(size = 22, hjust = 0.5),
          plot.title = element_text(size = 25, hjust = 0.5),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text = element_text(size = 15),
          # panel.background = element_rect(fill = "#f5f5f5", color = NA),
          # plot.background = element_rect(fill = "#f5f5f5", color = NA),
          # panel.grid.minor = element_line(colour = "darkgray"),
          # panel.grid.major = element_line(colour = "darkgray"),
          legend.text = element_text(size = 16),
          legend.title = element_blank()
        ) +
        scale_y_continuous(name = "Risk reduction")
    }
    else {
      ggplot(mapping = aes(x, 100*y)) +
        geom_point(size = 3, color = mode_color()) +
        geom_line(color = mode_color()) +
        # theme_minimal() +
        plot_theme() +
        labs(
          title = "Risk reduction",
          subtitle = subtitle,
          x = x_lab,
          y = "Risk reduction"
        ) +
        theme(
          plot.subtitle = element_text(size = 22, hjust = 0.5),
          plot.title = element_text(size = 25, hjust = 0.5),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text = element_text(size = 15),
          # panel.background = element_rect(fill = "#f5f5f5", color = NA),
          # plot.background = element_rect(fill = "#f5f5f5", color = NA),
          # panel.grid.minor = element_line(colour = "darkgray"),
          # panel.grid.major = element_line(colour = "darkgray")
        ) +
        scale_y_continuous(name = "Relative risk reduction (%)",
                           sec.axis = sec_axis(~.*ifelse(selected_x == "K", x, input$K/100), name = "Absolute risk reduction (%)")) 
    }
  }, message = "Calcuating...")})
  
  output$summary <- renderPrint({withProgress({
    if (!accepeted) showModal(modal)
    # print(input$K2)
    # cat("<div class = \"alert alert-info\">")
    
    # if (input$r2 >= input$h2) {
    #   updateNumericInput(session, "input_r2", value = input$h2)
    #   updateSliderInput(session, "r2", value = input$h2)
    # }
    
    # Possible states
    family_history <- input$siblings_n > 0 | input$parents_n > 0
    prs_condition <- input$type2 == "Conditional"
    exclude_strategy <- input$lowestexclude2 != "Lowest"
    
    binomial_model <- input$det_random2 == "Binomial"
    poisson_model <- input$det_random2 == "Poisson"
    
    if(binomial_model) {
      updateSliderInput(session, "N2", max = 20,
                        label = "Number of euploid embryos")
    }
    else if(poisson_model) {
      updateSliderInput(session, "N2", max = 20,
                        label = "Expected number of euploid embryos")
    }
    else {
      updateSliderInput(session, "N2", max = 10,
                        label = "Number of live births")
    }
    
    if(exclude_strategy) {
      updateSliderInput(session, "N2", max = 10,
                        label = "Number of live births")
    }
    
    if (!family_history) {
      # Either "normal" or "conditional"
      if(prs_condition) {
        if (exclude_strategy) {
          temp <-
            risk_reduction_exclude_conditional(
              input$r2,
              K = input$K2/100,
              q = input$q2,
              n = input$N2,
              qf = 1-input$qf2,
              qm = 1-input$qm2,
              relative = T)
        }
        else {
          temp <-
            risk_reduction_lowest_conditional(
              input$r2,
              K = input$K2/100,
              n = input$N2,
              qf = 1-input$qf2,
              qm = 1-input$qm2,
              relative = T,
              parental_avg_given = F)
          if (binomial_model) {
            temp <-
              risk_reduction_lowest_conditional_bin(
                input$r2,
                K = input$K2/100,
                n = input$N2,
                qf = 1-input$qf2,
                qm = 1-input$qm2,
                relative = T,
                parental_avg_given = F, input$p_lb2)
          }
          else if (poisson_model) {
            temp <-
              risk_reduction_lowest_conditional_pois(
                input$r2,
                K = input$K2/100,
                lambda = input$N2 * input$p_lb2,
                qf = 1-input$qf2,
                qm = 1-input$qm2,
                relative = T,
                parental_avg_given = F)
          }
        }
        print_result_conditional(temp)
      }
      else {
        if(exclude_strategy) {
          temp <- risk_reduction_exclude(input$r2,
                                         K = input$K2/100,
                                         q = input$q2,
                                         n = input$N2)
        }
        else {
          temp <- risk_reduction_lowest(input$r2, K = input$K2/100, 
                                        n = input$N2)
          if (binomial_model) {
            temp <- risk_reduction_lowest_bin(input$r2, K = input$K2/100, 
                                          n = input$N2, input$p_lb2)
          }
          else if (poisson_model) {
            temp <- risk_reduction_lowest_pois(input$r2, K = input$K2/100, 
                                              input$N2*input$p_lb2)
          }
        }
        print_result(temp, input$K2)
      }
    }
    else {
      # Condition on family history
      set.seed(1)
      if(prs_condition) {
        if (exclude_strategy) {
          temp <-
            risk_parents_offspring_generic(
              200000,
              input$N2,
              input$r2,
              input$h2,
              input$K2/100,
              input$sick_parents,
              input$parents_n-input$sick_parents,
              input$sick_siblings,
              input$siblings_n - input$sick_siblings, 
              1-input$qf2,
              1-input$qm2, 
              "exclude_percentile",
              input$q2)
          # temp <-
          #   risk_parents_offspring_generic(
          #     200000,
          #     input$N2,
          #     input$r2,
          #     input$h2,
          #     input$K2/100,
          #     input$sick_parents,
          #     input$parents_n-input$sick_parents,
          #     input$sick_siblings,
          #     input$siblings_n - input$sick_siblings, 
          #     1-input$qf2,
          #     1-input$qm2, 
          #     "exclude_percentile",
          #     input$q2, random_strategy = ifelse(binomial_model, "Binomial",
          #                      ifelse(poisson_model, "Poisson", "Fixed")),
          #     p = input$p_lb2)
        }
        else {
          temp <-
            risk_parents_offspring_generic(
              200000,
              input$N2,
              input$r2,
              input$h2,
              input$K2/100,
              input$sick_parents,
              input$parents_n-input$sick_parents,
              input$sick_siblings,
              input$siblings_n - input$sick_siblings,
              1-input$qf2,
              1-input$qm2, random_strategy = ifelse(binomial_model, "Binomial",
                                  ifelse(poisson_model, "Poisson", "Fixed")),
              p = input$p_lb2)
        }
      }
      else {
        if (exclude_strategy) {
          temp <-
            risk_parents_offspring_generic(
              200000,
              input$N2,
              input$r2,
              input$h2,
              input$K2/100,
              input$sick_parents,
              input$parents_n-input$sick_parents,
              input$sick_siblings,
              input$siblings_n - input$sick_siblings, 
              selection_strategy = "exclude_percentile",
              exclusion_q = input$q2)
          
          # temp <-
          #   risk_parents_offspring_generic(
          #     200000,
          #     input$N2,
          #     input$r2,
          #     input$h2,
          #     input$K2/100,
          #     input$sick_parents,
          #     input$parents_n-input$sick_parents,
          #     input$sick_siblings,
          #     input$siblings_n - input$sick_siblings, 
          #     selection_strategy = "exclude_percentile",
          #     exclusion_q = input$q2, random_strategy = ifelse(binomial_model, "Binomial",
          #                                                      ifelse(poisson_model, "Poisson", "Fixed")),
          #     p = input$p_lb2)
        }
        else {
          temp <-
            risk_parents_offspring_generic(
              200000,
              input$N2,
              input$r2,
              input$h2,
              input$K2/100,
              input$sick_parents,
              input$parents_n-input$sick_parents,
              input$sick_siblings,
              input$siblings_n - input$sick_siblings,
              random_strategy = ifelse(binomial_model, "Binomial",
                                         ifelse(poisson_model, "Poisson", "Fixed")),
              p = input$p_lb2)
        }
      }
      print_family_history(temp)
    }
    
    # cat("</div>")
    }, message = "Calcuating...")})
  
  output$two_traits <- renderPrint({
    cat("<div class = \"alert alert-info\">")
    temp <-
      simulate_lowest_risk_two_traits(input$r2_1, input$r2_2, input$rho, input$K_1, input$K_2/100, input$N_2, input$samples_2)
    cat(sprintf("<p><u>Relative risk reduction for disease 1</u>: <strong>%.4f</strong>\n</p>", temp[1]))
    cat(sprintf("<p><u>Absolute risk reduction for disease 1</u>: <strong>%.4f</strong>\n</p>", temp[2]))
    cat(sprintf("<p><u>Couples needed to screen for disease 1</u>: <strong>%.0f</strong>\n</p>", ceiling(1/temp[2])))
    
    cat(sprintf("<p><u>Relative risk reduction for disease 2</u>: <strong>%.4f</strong>\n</p>", ifelse(input$rho == 0, 0, temp[3])))
    cat(sprintf("<p><u>Absolute risk reduction for disease 2</u>: <strong>%.4f</strong>\n</p>", ifelse(input$rho == 0, 0, temp[4])))
    cat(sprintf("<p><u>Couples needed to screen for disease 2</u>: <strong>%.0f</strong>\n</p>", ifelse(input$rho == 0, 0, ceiling(1/temp[4]))))
    
    # cat(sprintf("<p style = \"color:red\">Based on %d simulations.</p>", input$samples_2))
    cat(sprintf("<p style = \"color:red\">Based on simulation.</p>", input$samples_2))
    cat("</div>")
  })
  
  output$logo_ui <- renderUI({
    logo_src <- if (input$dark_mode == "dark") {
      "Images/logo_dark.webp"  # Dark mode logo
    } else {
      "Images/logo.webp"        # Light mode logo
    }
    
    tags$a("", tags$img(src = logo_src, width = "200px", height = "80px"))
  })
}

# shinyApp(ui = ui, server = server)
# runApp(launch.browser = F)

