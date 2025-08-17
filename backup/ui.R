# TODO:
# Remove random number from calculator (include only in plots) - V
# Add Poisson - V
# Remove threshold and such plot - V
# Change names of things (in the conditional, in the "number of embryos" - number of births, number of euploids and such)

# Poisson/binomial - doesn't work for exclude without history

# 1. Add random number of embryos for all (seems fine)
# 2. Add conditional on parents' PRS + disease status - maybe works for only siblings?
# 3. Correct family history (when conditioning on one parent+sibling, and when you condition on PRS)
# 4. Switch the loading to https://shiny.posit.co/r/articles/build/progress/
# 5. Switch the rest of the UI to bslib as in the new calc (https://shiny.posit.co/r/articles/build/layout-guide/)
# 6. Have the plot and calculator modes synchronized - maybe not
# 7. Cookies to remember things between runs?
# 8. Cache for faster computations?
# 9. Minify things? (js, css)
# 10. Add dark mode?

library(shiny)
library(shinyWidgets)
library(bslib)
# library(thematic)

# thematic_shiny(font = 'auto')

slider_and_numeric <- function(id, label, min, max, step, value, helptext = "",
                               placement = "bottom",
                               post = NULL) {
  if(length(step) > 1) {
    div(id = id, 
        withMathJax(),
        splitLayout(cellWidths = c("80%", "20%"),
                    sliderTextInput(
                      inputId = id,
                      label = label,
                      choices = step,
                      grid = F, force_edges = T,
                      post = post),
                    helpPopup(NULL, helptext, placement = placement, c("hover"))))
  }
  else {
    div(id = id, 
    splitLayout(cellWidths = c("80%", "20%"),
                sliderInput(
                  inputId = id,
                  label = label,
                  min = min,
                  max = max,
                  step = step,
                  value = value,
                  post = post
                ), 
                helpPopup(NULL, helptext, placement = placement, c("hover"))))
  }
}


helpPopup <- function(title, content,
                      placement=c('right', 'top', 'left', 'bottom'),
                      trigger=c('click', 'hover', 'focus', 'manual')) {
  tooltip(bsicons::bs_icon("question-circle"), content, id="tip", placement = placement)
}

disclamir_and_date_text <- HTML("<p align=\"center\"><b><font color = \"red\">The application is intended for research purposes only and is not intended to guide clinical decision making</font></b><br>",
                                paste("Last update", max(format(max(file.info("ui.R")$mtime,
                                                                    file.info("server.R")$mtime,
                                                                    file.info("EmbryoSelection.R")$mtime), "%d-%m-%Y"))),
                                "</p>")

# Another option?
about_panel <- page_fillable(h1("About", align = "center"),
                                  HTML("<p>This application allows users to predict the expected risk reduction when selecting an IVF embryo for transfer based on polygenic risk scores (PRS) for a single disease.</p>"),
                                  p("We provide estimates under two selection strategies."),
                                  HTML("<ol>
                                   <li>Lowest risk prioritization: Selecting the embryo with the lowest PRS among all available embryos.</li>
                                   <li>High-risk exclusion: Excluding embryos with a PRS above a “high-risk” cutoff, and then selecting an embryo at random among the remaining embryos. In the case all embryos are high-risk, a random embryo is selected.</li></ol>"),
                                  h2("Plot mode"),
                                  p("Use the Plot mode to generate graphs of the relative and absolute risk reductions vs the parameters of the problem for each of the two selection strategies."),
                                  p("The parameters are:"),
                                  HTML("<ol>
                                   <li>PRS accuracy ($R^2$): The proportion of variance in liability to the disease explained by the PRS, which is a measure of the accuracy of the score. $R^2$ is currently between 5-10% for most common polygenic diseases.</li>
                                   <li>The disease prevalence: The proportion of individuals in the (adult) population affected by the disease.</li>
                                   <li>The number of embryos: The number of viable embryos from which a single embryo is selected for transfer. Important note: the model assumes that each embryo transferred will be born. This parameter should therefore correspond to the number of live births expected from the given cycle.</li>
                                   <li>Quantile from which to exclude: In the “high-risk exclusion” strategy, this is the cutoff that defines embryos as high-risk. Embryos with PRS above that quantile are excluded. For example, if the parameter equals 10%, all embryos with PRS at the top 10% of the distribution of PRS in the population will be excluded.</li>
                                   </ol>"),
                                  h2("Calculator mode"),
                                  p("Use the calculator mode to present the baseline risk and the risk reduction given a set of parameters. The calculator outputs the risk with and without embryo selection, the relative and the absolute risk reduction, and the number of couples needed to screen their embryos to prevent a single future case. The calculator also allows users to compute risk estimates when conditioning on the parental PRS percentiles or disease status."),
                                  p("When conditioning on the parental PRS percentiles, the following parameters are required:"),
                                  HTML("<ol>
                                   <li>Father's PRS percentile: For example, if this parameter equal 0.05, the PRS of the father is at the top 5% of the distribution of PRS in the population.</li>
                                   <li>Mother's PRS percentile: An analogous parameter for the mother.</li>
                                   </ol>"),
                                  p("When conditioning on the parental disease status, the following parameters are required:"),
                                  # HTML("<ol>
                                  #      <li>Father has the disease: yes/no</li>
                                  #      <li>Mother has the disease: yes/no</li>
                                  #      <li>$h^2$: The heritability of the disease, i.e., the proportion of variance in the liability of the disease explained by additive genetic factors.</li>
                                  #      <li>Number of Monte-Carlo draws: when conditioning on the parental disease status, the risk is estimated using a Monte-Carlo method. With more draws, the risk estimate becomes more accurate but is slower to compute.</li>
                                  #      </ol>"),
                                  HTML("<ol>
                                   <li>Father has the disease: yes/no</li>
                                   <li>Mother has the disease: yes/no</li>
                                   <li>$h^2$: The heritability of the disease, i.e., the proportion of variance in the liability of the disease explained by additive genetic factors.</li>
                                   </ol>"),
                                  # h2("Two diseases"),
                                  # p("Two diseases is a calculator of the risk of two diseases, given a certain correlation and when selecting the embryo with the lowest PRS for disease 1. It is based on a simulation as explained in the paper."),
                                  # p("The only new variable here is $\\rho$, the genetic correlation between the two diseases."),
                                  h2("Reference"),
                                  HTML("<ol>
                                   <li> <cite>Lencz, T., Backenroth, D., Granot-Hershkovitz, E., Green, A., Gettler, K., Cho, J. H., Weissbrod, O., Zuk, O., & Carmi, S. (2021). Utility of polygenic embryo screening for disease depends on the selection strategy. eLife, 10. <a href=\"https://doi.org/10.7554/elife.64716\" target=\"_blank\" rel=\"noopener noreferrer\">https://doi.org/10.7554/elife.64716</a> </cite> </li>
                                   <li> <cite>Capalbo A, de Wert G, Mertes H, Klausner L, Coonen E, Spinella F, Van de Velde H, Viville S, Sermon K, Vermeulen N, Lencz T, Carmi S. Screening embryos for polygenic disease risk: a review of epidemiological, clinical, and ethical considerations. Hum Reprod Update. 2024 Oct 1; <a href=\"https://doi.org/10.1093/humupd/dmae012\" target=\"_blank\" rel=\"noopener noreferrer\">https://doi.org/10.1093/humupd/dmae012</a> </cite> </li>
                                   <li> <cite>Klausner, L., Lencz, T. & Carmi, S. (2025). PEStimate: Predicting offspring disease risk after Polygenic Embryo Screening.</cite> </li>
                                   </ol>"),
                                  h2("Contact"),
                                  p("Please contact us if you find an error or have any suggestion."),
                                  HTML("<p>Shai Carmi, <a href=\"mailto: shai.carmi@huji.ac.il\">shai.carmi@huji.ac.il</a></p>"),
                                  HTML("<p>Liraz Klausner, <a href=\"mailto: liraz.klausner@mail.huji.ac.il\">liraz.klausner@mail.huji.ac.il</a></p>"),
                                  p("Braun School of Public Health, The Hebrew University of Jerusalem"),
                                  disclamir_and_date_text)

# plot_panel <- page_fillable(
#   layout_columns(card(
#     conditionalPanel("input.lowestexclude == \"Lowest\"",
#                      selectInput("det_random", label = "Is the number of embryos fixed?",
#                                  choices = c("Fixed live births", 
#                                              "Binomial",
#                                              "Poisson"),
#                                  width = "50%"),
#                      conditionalPanel("input.det_random != \"Fixed live births\"", 
#                                       slider_and_numeric("p_lb", "Probability of live birth:",
#                                                          0.01, 0.99, 0.01, value=0.3, 
#                                                          helptext = "Out of the available euploid embryos, what proportion of them would result in a live birth?"))),
#     radioButtons("x_var", "Variable for x axis", choiceNames = c("R-squared", "Disease prevalence", "Number of embryos"), 
#                  choiceValues = c("r2", "Disease prevalence", "Number of embryos"),
#                  selected = "r2",
#                  inline = T),
#     radioButtons(inputId = "lowestexclude",
#                  label = "Choose lowest risk embryo or exclude high risk embroys",
#                  choices = c("Lowest", "Exclude"), inline = T
#     )
#   ),
#   card(
#     conditionalPanel("input.x_var != \"r2\"",
#                      slider_and_numeric("r", "PRS accuracy ($R^2$):", 0.01, 1, NULL, 0.05, 
#                                         "The proportion of the variance in the liability of the disease explained by the polygenic risk score. It is a measure of the accuracy of the score. Typically in the range 0.05-0.1.")),
#     conditionalPanel("input.x_var != \"Disease prevalence\"",
#                      slider_and_numeric("K", "Disease prevalence:", 0.001*100, 0.3*100, 
#                        100*sort(unique(c(seq(0.001, 0.3, 0.001), 
#                                          round(exp(seq(log(0.001), log(0.3), length = 500)), digits = 4)))), 100*0.001, 
#                        "What percent of the population have the disease?",
#                        post = "%")),
#     conditionalPanel("input.x_var != \"Number of embryos\"",
#                      slider_and_numeric("N", "Number of embryos:", 2, 10, 1, 5, "The number of embryos available for selection.")),
#     conditionalPanel(
#       condition = "input.lowestexclude == 'Exclude' & input.x_var != 'Percentile'",
#       slider_and_numeric("q", "Percentile from which to exclude embryos:", 0.01, 0.99, 0.01, 0.3, paste("Embryos with PRS above that percentile are excluded. For example, if the parameter equals 0.1, all embryos with PRS at the top 10% of the distribution of PRS in the population will be excluded. If no embryo is below the threshold we select one randomly.")))
#   )),
#   card(plotOutput(outputId = "distPlot")),
#   disclamir_and_date_text
# )

plot_panel <- page_fillable(
  layout_columns(col_widths = c(4, 8), fill = F,
                 div(
  card(
    conditionalPanel("input.lowestexclude == \"Lowest\"",
                     selectInput("det_random", label = "Is the number of embryos fixed?",
                                 choices = c("Fixed live births", 
                                             "Binomial",
                                             "Poisson"),
                                 width = "50%"),
                     conditionalPanel("input.det_random != \"Fixed live births\"", 
                                      slider_and_numeric("p_lb", "Probability of live birth:",
                                                         0.01, 0.99, 0.01, value=0.3, 
                                                         helptext = "Live birth rate per euploid embryo transfer"))),
    radioButtons("x_var", "Variable for x axis", choiceNames = c("R-squared", "Disease prevalence", "Number of embryos"), 
                 choiceValues = c("r2", "Disease prevalence", "Number of embryos"),
                 selected = "r2",
                 inline = T),
    radioButtons(inputId = "lowestexclude",
                 label = "Choose lowest risk embryo or exclude high risk embroys",
                 choices = c("Lowest", "Exclude"), inline = T
    )),
  card(
    conditionalPanel("input.x_var != \"r2\"",
                     slider_and_numeric("r", "PRS accuracy ($R^2$):", 0.01, 1, NULL, 0.05, 
                                        "The proportion of the variance in the liability of the disease explained by the polygenic risk score. It is a measure of the accuracy of the score. Typically in the range 0.05-0.1.")),
    conditionalPanel("input.x_var != \"Disease prevalence\"",
                     slider_and_numeric("K", "Disease prevalence:", 0.001*100, 0.3*100, 
                                        100*sort(unique(c(seq(0.001, 0.3, 0.001), 
                                                          round(exp(seq(log(0.001), log(0.3), length = 500)), digits = 4)))), 100*0.001, 
                                        "Fraction of the population with the disease",
                                        post = "%")),
    conditionalPanel("input.x_var != \"Number of embryos\"",
                     slider_and_numeric("N", "Number of embryos:", 2, 10, 1, 5, "The number of embryos/births avaliable for selection.")),
    conditionalPanel(
      condition = "input.lowestexclude == 'Exclude' & input.x_var != 'Percentile'",
      slider_and_numeric("q", "Percentile from which to exclude embryos:", 0.01, 0.99, 0.01, 0.3, paste("Embryos with PRS above that percentile are excluded. For example, if the parameter equals 0.1, all embryos with PRS at the top 10% of the distribution of PRS in the population will be excluded. If no embryo is avaliable, select one at random.")))
  )),
  card(plotOutput(outputId = "distPlot"))),
  disclamir_and_date_text
)

calc_panel <- page_fillable(
  tags$head(tags$script(("
    Number.prototype.mapLog = function (min, max) {
      const mapped = (this - min) * (Math.log(max) - Math.log(min)) / (max - min) + Math.log(min);
      return Math.exp(mapped);
    }
    $(function() {
    setTimeout(function(){
    
    # $('#K2').data('ionRangeSlider').update({
    #        'prettify': function (n) {
    #        return (n.mapLog(this.min, this.max).toLocaleString('en-US'));
    #        }})

    $('#K_1').data('ionRangeSlider').update({
           'prettify': function (n) {
           return (n.mapLog(this.min, this.max).toLocaleString('en-US'));
           }})
    $('#K_2').data('ionRangeSlider').update({
           'prettify': function (n) {
           return (n.mapLog(this.min, this.max).toLocaleString('en-US'));
           }})
                           }, 2)})"))),
  layout_columns(col_widths = c(4, 8), fill=F, 
                 # layout_column_wrap(heights_equal = "row", width=1,
                 div(
                   card(card_header("Embryo parameters"),
                        conditionalPanel("input.lowestexclude2 == \"Lowest\"",
                                         selectInput("det_random2", label = "Is the number of embryos fixed?",
                                                     choices = c("Fixed live births", 
                                                                 "Binomial",
                                                                 "Poisson"),
                                                     width = "50%"),
                                         conditionalPanel("input.det_random2 != \"Fixed live births\"", 
                                                          slider_and_numeric("p_lb2", "Probability of live birth:",
                                                                             0.01, 0.99, 0.01, value=0.3,
                                                                             helptext = "Live birth rate per euploid embryo transfer"))),
                        slider_and_numeric("N2", "Number of live births:", 2, 10, 1, 5, "The number of embryos/births avaliable for selection."),
                        radioButtons(
                          inputId = "lowestexclude2",
                          label = "Choose lowest risk embryo or exclude high risk embroys",
                          choices = c("Lowest", "Exclude"), inline = T)),
                   card(card_header("Disease parameters"),
                        # slider_and_numeric("K2", "Disease prevalence:", 0.01, 1, 0.01, 0.5, NULL),
                        slider_and_numeric("K2", "Disease prevalence:", 100*0.001, 100*0.3, 
                                           100*sort(unique(c(seq(0.001, 0.3, 0.001), 
                                                             round(exp(seq(log(0.001), log(0.3), length = 500)), digits = 4)))), 100*0.001, "Fraction of the population with the disease",
                                           post = "%"),
                        slider_and_numeric("r2", "PRS accuracy ($R^2$):", 0.01, 0.99, NULL, 0.05, "The proportion of the variance in the liability of the disease explained by the polygenic risk score. It is a measure of the accuracy of the score. Typically in the range 0.05-0.1. Must be smaller than $h^2$ when conditioning on family status."),
                        slider_and_numeric("q2", "Percentile from which to exclude embryos:", 0.01, 0.99, 0.01, 0.3, paste("Embryos with PRS above that percentile are excluded. For example, if the parameter equals 0.1, all embryos with PRS at the top 10% of the distribution of the PRS in the population will be excluded. If no embryo is avaliable, select one at random.")),
                        slider_and_numeric("h2", "$h^2:$", 0.01, 0.99, 0.01, 0.4, "The heritability of the disease. Only relevant when conditioning on the family disease status."))),
                 card(card_header("Family information", style = "text-align: center;"),
                      layout_column_wrap(width=1/2,
                                         card(card_header("Family disease status"),
                                           # conditionalPanel(
                                           # condition = "input.type2 == 'Conditional'",
                                           # ),
                                           # conditionalPanel(
                                           # condition = "input.type2 == 'Family History'",
                                           sliderInput(inputId ="parents_n",
                                                       label = "Number of parents with known disease status:",
                                                       min=0,
                                                       max=2,
                                                       value=0,
                                                       step=1),
                                           sliderInput(inputId = "sick_parents",
                                                       label = "Number of affected parents:",
                                                       min = 0,
                                                       max = 2,
                                                       value = 0,
                                                       step = 1),
                                           sliderInput("siblings_n",
                                                       label = "Number of sibling with known disease status:",
                                                       min = 0,
                                                       max = 20,
                                                       value = 0,
                                                       step = 1),
                                           sliderInput("sick_siblings",
                                                       label = "Number of affected siblings:",
                                                       min = 0,
                                                       max = 20,
                                                       value = 0,
                                                       step = 1)),
                                         div(
                                             card(card_header("Parental polygenic risk score"),
                                               radioButtons(
                                                 inputId = "type2",
                                                 label = "Condition on the parentel polygenic risk score?",
                                                 # choices = c("Risk reduction", "Conditional", "Family History"), inline = T
                                                 # choiceValues = c("Risk reduction", "Conditional", "Family History"),
                                                 # choiceNames = c("No conditioning", "Conditional on the parents' polygenic risk score", "Conditional on family disease status"),
                                                 choiceValues = c("Risk reduction", "Conditional"),
                                                 choiceNames = c("No", "Yes"), inline = T,
                                               ),
                                               slider_and_numeric("qf2", "Father's polygenic risk score percentile:", 0.01, 0.99, 0.01, 0.5, paste("For example, if this parameter equal 0.05, the PRS of the father is at the top 5% of the distribution of the PRS in the population.")),
                                               slider_and_numeric("qm2", "Mother's polygenic risk score percentile:", 0.01, 0.99, 0.01, 0.5, paste("For example, if this parameter equal 0.05, the PRS of the mother is at the top 5% of the distribution of the PRS in the population."))), 
                                             card(card_header("Results"),
                                                  class = "alert alert-info", htmlOutput("summary")))))
  ),
  disclamir_and_date_text
)

calc_two_traits <- div(class = "well", fluidRow(column(4,
                                                       slider_and_numeric("N_2", "Number of embryos:", 2, 10, 1, 5, "The number of embryos available for selection."),
                                                       slider_and_numeric("rho", '$\\rho$, the genetic correlation between the diseases:', -0.99, 0.99, 0.01, 0, "The genetic correlation between the two diseases."),
                                                       slider_and_numeric("samples_2", "Number of monte carlo draws:", 100000, 500000, 1000, 100000, "The number of simulations. Higher number will give a more accurate estimate, but might take longer to run.")),
                                                column(4, 
                                                       slider_and_numeric("r2_1", "PRS accuracy ($R^2 ~ \\text{disease 1}$):", 0.01, 1, 0.001, 0.05, "The proportion of the variance in the liability of the first disease explained by the polygenic risk score. It is a measure of the accuracy of the score. Typically in the range 0.05-0.1."),
                                                       slider_and_numeric("r2_2", "PRS accuracy ($R^2 ~ \\text{disease 2}$):", 0.01, 1, 0.001, 0.05, "The proportion of the variance in the liability of the second disease explained by the polygenic risk score. It is a measure of the accuracy of the score. Typically in the range 0.05-0.1."),
                                                       fluidRow(column(8, offset = 2, htmlOutput("two_traits"), align = "center"))),
                                                column(4, 
                                                       slider_and_numeric("K_1", "Prevalence of disease 1:", 0.001, 0.3, unique(round(exp(seq(log(0.001), log(0.3), length = 500)), digits = 4)), 0.001, "How prevalent is the first disease in the population? 0.01 means that 1% of the population have the disease, and 0.2 means that 20% of the population have the disease."),
                                                       slider_and_numeric("K_2", "Prevalence of disease 2:", 0.001, 0.3, unique(round(exp(seq(log(0.001), log(0.3), length = 500)), digits = 4)), 0.001, "How prevalent is the second disease in the population? 0.01 means that 1% of the population have the disease, and 0.2 means that 20% of the population have the disease."))),
                       disclamir_and_date_text)



# addResourcePath(prefix = "imgResources", directoryPath = "myimages")
addResourcePath(prefix = "Images", directoryPath = "Images/")

# ui <- fluidPage(
ui <- page_navbar(
  # header = tags$head(
    # tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-AMS_HTML")
  # ),
  # tags$head(HTML("<script type=\"text/x-mathjax-config\">
  # MathJax.Hub.Config({
  #   tex2jax: {
  #     inlineMath: [ ['$','$'] ],
  #     processEscapes: true
  #   }
  # });
  # </script>"),
  #           # Thanks to https://stackoverflow.com/questions/17325521/r-shiny-display-loading-message-while-function-is-running
  #           tags$style(type="text/css", "
  #          #loadmessage {
  #            position: fixed;
  #            top: 0px;
  #            left: 0px;
  #            width: 100%;
  #            padding: 5px 0px 5px 0px;
  #            text-align: center;
  #            font-weight: bold;
  #            font-size: 100%;
  #            color: #000000;
  #            background-color: #CCFF66;
  #            z-index: 105;
  #          }")),
  
  # shinyjs::useShinyjs(),
  header=withMathJax(tags$head(HTML("<script type='text/x-mathjax-config' >
    MathJax.Hub.Config({
    tex2jax: {inlineMath: [['$','$']]}
    });
    </script>")),
                     shinyjs::useShinyjs()),
  # tags$div("Here is an equation: $E = mc^2$"),
  # input_dark_mode(id = "mode"),
  # tabsetPanel(
  #   type = "pills",
  #   nav_item(input_dark_mode()),
  #   tabPanel("About", about_panel),
  #   tabPanel("Plot", plot_panel),
  #   tabPanel("Calculator", calc_panel),
  #   # tabPanel("Two diseases", calc_two_traits)
  # )
  # title="PEStimate",
  title=tags$a("", uiOutput("logo_ui")),
  nav_panel("About", about_panel),
  nav_panel("Plot", plot_panel),
  nav_panel("Calculator", calc_panel),
  nav_spacer(),
  nav_item(input_dark_mode(id = "dark_mode")),
)
