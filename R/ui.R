# TODO (?):
# Have the plot and calculator modes synchronized - maybe not
# Cookies to remember things between runs?
# Cache results?
# Minify things?

slider_and_numeric <- function(id, label, min, max, step, value, helptext = "",
                               placement = "bottom",
                               post = NULL) {
  if(length(step) > 1) {
    sliderTextInput(inputId = id,
                    # label = div(class = "position-relative d-flex align-items-center",
                    label = span(
                      label,
                      helpPopup(NULL, helptext, placement = placement, c("hover"))),
                    choices = step,
                    grid = F, force_edges = T,
                    post = post)
    # div(id = id, 
    #     withMathJax(),
    #     splitLayout(cellWidths = c("80%", "20%"),
    #                 sliderTextInput(
    #                   inputId = id,
    #                   label = label,
    #                   choices = step,
    #                   grid = F, force_edges = T,
    #                   post = post),
    #                 helpPopup(NULL, helptext, placement = placement, c("hover"))))
  }
  else {
    sliderInput(
      inputId = id,
      label = span(label,
                   helpPopup(NULL, helptext, placement = placement, c("hover"))),
      min = min,
      max = max,
      step = step,
      value = value,
      post = post)
    
    # div(id = id, 
    # splitLayout(cellWidths = c("80%", "20%"),
    #             sliderInput(
    #               inputId = id,
    #               label = label,
    #               min = min,
    #               max = max,
    #               step = step,
    #               value = value,
    #               post = post
    #             ), 
    #             helpPopup(NULL, helptext, placement = placement, c("hover"))))
  }
}


helpPopup <- function(title, content,
                      placement=c('right', 'top', 'left', 'bottom'),
                      trigger=c('click', 'hover', 'focus', 'manual')) {
  # tooltip(bsicons::bs_icon("question-circle"), content, id="tip", placement = placement)
  tooltip(HTML('<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" fill="currentColor" class="bi bi-question-circle position-abolute end-0" viewBox="0 0 16 16">
  <path d="M8 15A7 7 0 1 1 8 1a7 7 0 0 1 0 14m0 1A8 8 0 1 0 8 0a8 8 0 0 0 0 16"/>
  <path d="M5.255 5.786a.237.237 0 0 0 .241.247h.825c.138 0 .248-.113.266-.25.09-.656.54-1.134 1.342-1.134.686 0 1.314.343 1.314 1.168 0 .635-.374.927-.965 1.371-.673.489-1.206 1.06-1.168 1.987l.003.217a.25.25 0 0 0 .25.246h.811a.25.25 0 0 0 .25-.25v-.105c0-.718.273-.927 1.01-1.486.609-.463 1.244-.977 1.244-2.056 0-1.511-1.276-2.241-2.673-2.241-1.267 0-2.655.59-2.75 2.286m1.557 5.763c0 .533.425.927 1.01.927.609 0 1.028-.394 1.028-.927 0-.552-.42-.94-1.029-.94-.584 0-1.009.388-1.009.94"/>
</svg>'), content, id="tip", placement = placement)
}

disclamir_and_date_text <- HTML("<p align=\"center\"><b><font color = \"red\">The application is intended for research purposes only and is not intended to guide clinical decision making</font></b><br>",
                                paste("Last update", max(format(max(file.info("R/ui.R")$mtime,
                                                                    file.info("R/server.R")$mtime,
                                                                    file.info("R/EmbryoSelection.R")$mtime), "%d-%m-%Y"))),
                                "</p>")


about_panel <- page_fillable(h1("About", align = "center"),
                             # HTML("<p>This application allows users to predict the expected risk reduction when selecting an IVF embryo for transfer based on polygenic risk scores (PRS) for a single disease.</p>"),
                             p("The app predicts the risk reduction when selecting an IVF embryo for transfer based on polygenic risk scores (PRS) for a single disease."),
                             HTML("<p>Use the <strong>Calculator tab</strong> to present various metrics for risk reduction given any set of model parameters.<br>
                               Use the <strong>Plot tab</strong> to generate graphs of the risk reductions vs some of the parameters.</p>"),
                             h2("Model parameters"),
                             h3("Selection strategy"),
                             p("Users should select one of the following embryo selection strategies."),
                             HTML("<ol>
                                   <li>Lowest risk prioritization. Select the embryo with the lowest PRS among all available embryos.</li>
                                   <li>High-risk exclusion. Exclude embryos with a PRS above a “high-risk” cutoff (to be specified), and then select an embryo at random from among the remaining embryos. In case all embryos are high-risk, select a random embryo.</li></ol>"),
                             h3("Embryo parameters"),
                             p("Users should select one of the following models and specify the associated parameters."),
                             HTML("<ol>
                                  <li>A fixed number of live births. There are n euploid embryos, and each will be born when transferred.</li>
                                  <li>A binomial number of live births: There are n euploid embryos, and each will be born with probability p.</li>
                                  <li>A Poisson number of live births: The number of euploid embryos has a Poisson distribution with mean n, and each will be born with probability p.</li></ol>"),
                             h3("Information on family members"),
                             p("Users can specify the disease status and genetic risk of family members."),
                             HTML("<ol>
                                   <li>Whether the <b>grandparents</b> from either side are <b>sick</b>, <b>healthy</b> or with an <b>unknown</b> disease status.</li>
                                   <li>Whether the <b>parents</b> from either side are <b>sick</b>, <b>healthy</b> or with an <b>unknown</b> disease status.</li>
                                   <li>The PRS percentile of each <b>parent</b>.</li>
                                   <li>Number of <b>siblings of each parent</b> with a known disease status.</li>
                                   <li>Number of <b>siblings of each parent</b> affected by the disease.</li>
                                   <li>Number of <b>siblings of the embryos</b> with a known disease status.</li>
                                   <li>Number of <b>siblings of the embryos</b> affected by the disease.</li>
                                   </ol>"),
                             h3("Disease parameters"),
                             p("Users should specify the following parameters."),
                             HTML("<ol>
                                   <li>PRS accuracy ($r^2$). The proportion of variance in liability to the disease explained by the PRS. This is a commonly used measure of PRS accuracy. $r^2$ is currently between 5-20% for most common polygenic diseases.</li>
                                   <li>$r^2$ adjustments. Users can adjust $r^2$ due to the target population.</li>
                                   <li>The disease prevalence. The proportion of individuals in the (adult) population affected by the disease.</li>
                                   <li>The heritability of the disease ($h^2$). The proportion of variance in the liability of the disease explained by additive genetic factors. (Only required when specifying the disease status of family members.)</li>
                                   </ol>"),
                             p("Users can also select from predefined diseaes, with estimates taken from the literature."),
                             # h3("$r^2$ adjustments"),
                             # p("Users can adjust the $r^2$ due to being from a different population, or due to $r^2$ being measured in a different scale then the liability scale."),
                             h3("Output"),
                             p("The app reports the following outputs."),
                             HTML("<ol>
                                   <li>The baseline risk (i.e., the risk of a randomly selected embryo).</li>
                                   <li>The risk of the embryo selected based on PRS.</li>
                                   <li>The probability of no risk reduction (one or no live births) when the number of live births is random.</li>
                                   <li>The relative risk reduction.</li>
                                   <li>The absolute risk reduction.</li>
                                   </ol>"),
                             h2("References"),
                             HTML("<ol>
                                   <li> <cite>Lencz, T., Backenroth, D., Granot-Hershkovitz, E., Green, A., Gettler, K., Cho, J. H., Weissbrod, O., Zuk, O., & Carmi, S. (2021). Utility of polygenic embryo screening for disease depends on the selection strategy. eLife. <a href=\"https://doi.org/10.7554/elife.64716\" target=\"_blank\" rel=\"noopener noreferrer\">https://doi.org/10.7554/elife.64716</a> </cite> </li>
                                   <li> <cite>Capalbo, A., de Wert, G., Mertes, H., Klausner, L., Coonen, E., Spinella, F., Van de Velde, H., Viville, S., Sermon, K., Vermeulen, N., Lencz, T. & Carmi, S. (2024) Screening embryos for polygenic disease risk: a review of epidemiological, clinical, and ethical considerations. Hum Reprod Update. <a href=\"https://doi.org/10.1093/humupd/dmae012\" target=\"_blank\" rel=\"noopener noreferrer\">https://doi.org/10.1093/humupd/dmae012</a> </cite> </li>
                                   <li> <cite>Klausner, L., Revital, A., Lencz, T. & Carmi, S. (2025). PEStimate: Predicting offspring disease risk after Polygenic Embryo Screening.</cite> </li>
                                   </ol>"),
                             h2("Contact"),
                             p("Please contact us if you find an error or have any suggestion."),
                             HTML("<p>Shai Carmi, <a href=\"mailto: shai.carmi@huji.ac.il\">shai.carmi@huji.ac.il</a></p>"),
                             HTML("<p>Liraz Klausner, <a href=\"mailto: liraz.klausner@mail.huji.ac.il\">liraz.klausner@mail.huji.ac.il</a></p>"),
                             p("Braun School of Public Health, The Hebrew University of Jerusalem"),
                             disclamir_and_date_text)

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
                                                                          0.01, 0.99, 0.01, value=0.5,
                                                                          helptext = "Live birth rate per euploid embryo transfer"))),
                     radioButtons("x_var", "Variable for x axis", choiceNames = c("r²", "Disease prevalence", "Number of embryos"), 
                                  choiceValues = c("r2", "Disease prevalence", "Number of embryos"),
                                  selected = "r2",
                                  inline = T),
                     radioButtons(inputId = "lowestexclude",
                                  label = "Choose lowest risk embryo or exclude high risk embroys",
                                  choices = c("Lowest", "Exclude"), inline = T
                     )),
                   card(
                     conditionalPanel("input.x_var != \"r2\"",
                                      slider_and_numeric("r", "PRS accuracy ($r^2$):", 0.01, 1, NULL, 0.05, 
                                                         "The proportion of the variance in the liability of the disease explained by the polygenic risk score. It is a measure of the accuracy of the score. Typically in the range 0.05-0.15.")),
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
                       slider_and_numeric("q", "Percentile from which to exclude embryos:", 0.01, 0.99, 0.01, 0.3, paste("Embryos with PRS above that percentile are excluded. For example, if the parameter equals 0.9, all embryos with PRS at the top 10% of the distribution of PRS in the population will be excluded. If no embryo is avaliable, select one at random.")))
                   )),
                 card(plotOutput(outputId = "distPlot"))),
  disclamir_and_date_text
)

embryo_parameters <- card(card_header("Embryo parameters", class = "bg-primary"),
                          conditionalPanel("input.lowestexclude2 == \"Lowest\"",
                                           selectInput("det_random2", label = "Is the number of embryos fixed?",
                                                       choices = c("Fixed live births", 
                                                                   "Binomial",
                                                                   "Poisson"),
                                                       width = "50%"),
                                           conditionalPanel("input.det_random2 != \"Fixed live births\"", 
                                                            slider_and_numeric("p_lb2", "Probability of live birth:",
                                                                               0.01, 0.99, 0.01, value=0.48,
                                                                               helptext = "Live birth rate per euploid embryo transfer"))),
                          slider_and_numeric("N2", "Number of live births:", 2, 10, 1, 5, "The number of embryos/births avaliable for selection."),
                          radioButtons(
                            inputId = "lowestexclude2",
                            label = "Choose lowest risk embryo or exclude high risk embroys",
                            choices = c("Lowest", "Exclude"), inline = T),
                          slider_and_numeric("q2", "Percentile from which to exclude embryos:", 0.01, 0.99, 0.01, 0.3, paste("Embryos with PRS above that percentile are excluded. For example, if the parameter equals 0.9, all embryos with PRS at the top 10% of the distribution of the PRS in the population will be excluded. If no embryo is avaliable, select one at random.")))

disease_parameters <- card(card_header("Disease parameters", class = "bg-primary"),
                           selectInput("disease_presets", "Disease Presets", 
                                       # choices = c("Custom", "Type 1 diabetes"),
                                       choices = "Custom",
                                       selected = "Custom"),
                           # slider_and_numeric("K2", "Disease prevalence:", 0.01, 1, 0.01, 0.5, NULL),
                           slider_and_numeric("K2", "Disease prevalence:", 100*0.001, 100*0.5, 
                                              100*sort(unique(c(seq(0.001, 0.5, 0.001), 
                                                                round(exp(seq(log(0.001), log(0.5), length = 500)), digits = 4)))), 100*0.001, "Fraction of the population with the disease",
                                              post = "%"),
                           slider_and_numeric("r2", "PRS accuracy ($r^2$):", 0.01, 0.99, NULL, 0.05, "The proportion of the variance in the liability of the disease explained by the polygenic risk score. It is a measure of the accuracy of the score. Typically in the range 0.05-0.15. Must be smaller than $h^2$ when conditioning on family disease status."),
                           selectizeInput("pop_adjust", span("Adjust $r^2$ for population?",
                                                          helpPopup(NULL, "Adjusts the $r^2$ based on estimates of partial correlation between UK biobank and other population", "bottom", "hover")),
                                       # tooltip(trigger = bsicons::bs_icon("question-circle-fill", size = "0.8em", style = "cursor: pointer; color: #007bc2;"),
                                       # "Adjusts the $r^2$ based on estimates of partial correlation between UK biobank and other population")),
                                       choices = c("Original"),
                                       options = list(dropdownParent = 'body')),
                           htmlOutput("r2_adjusted_output"),
                           slider_and_numeric("h2", "$h^2:$", 0.01, 0.99, 0.01, 0.4, "The heritability of the disease. Only relevant when conditioning on the family disease status."))

# r2_adj <- card(card_header("$r^2$ adjustments", class = "bg-primary"),
#                selectInput("pop_adjust", span("Adjust $r^2$ for population?",
#                                               helpPopup(NULL, "Adjusts the $r^2$ based on estimates of partial correlation between UK biobank and other population", "bottom", "hover")),
#                            # tooltip(trigger = bsicons::bs_icon("question-circle-fill", size = "0.8em", style = "cursor: pointer; color: #007bc2;"),
#                            # "Adjusts the $r^2$ based on estimates of partial correlation between UK biobank and other population")),
#                            choices = c("Original")),
#                
#                # splitLayout(cellWidths = c("80%", "20%"),
#                #             selectInput("pop_adjust", "Adjust $r^2$ for population?",
#                #                         choices = c("Original")), 
#                #             helpPopup(NULL, "Adjusts the $r^2$ based on estimates of partial correlation between UK biobank and other population", "bottom", "hover")),
#                slider_and_numeric("r2_adjust", "Adjust $r^2$ for lower liability variance explained?", 0.01, 1, 0.01, 1, helptext = "If the reported $r^2$ is given for the population, and you want to adjust it to a lower liability $r^2$."),
#                htmlOutput("r2_adjusted_output"))
#                # HTML("<p>Note: r² adjustment due to population are based on estimated correlation between UK and other populations.</p>"))

side <- function(num) {
  card(card_header(sprintf("Parent %d side", num), class = "bg-secondary text-white"),
       card_body(
         h5("Grandparents"),
         div(class = "d-flex gap-2",
             selectInput(sprintf("gp%da_status", num), "Grandparent 1", 
                         choices = c("Unknown"=0, "Sick"=1, "Healthy"=-1), width = "50%"),
             selectInput(sprintf("gp%db_status", num), "Grandparent 2", 
                         choices = c("Unknown"=0, "Sick"=1, "Healthy"=-1), width = "50%")
         ),
         h5(sprintf("Parent %d", num)),
         selectInput(sprintf("p%d_status", num), "Status", choices = c("Unknown"=0, "Sick"=1, "Healthy"=-1)),
         # Collapsible PRS for Parent 1
         # accordion(id="acc1",
         #   accordion_panel(
         #     "Known PRS (Optional)",
         #     value = "prs_panel",
         #     slider_and_numeric("qf2", "Parent 1 polygenic risk score percentile:", 0.01, 0.99, 0.01, 0.5, paste("For example, if this parameter equal 0.95, the PRS of the father is at the top 5% of the distribution of the PRS in the population.")),
         #   ),
         #   open = FALSE
         # ),
         radioButtons(
           inputId = "type2",
           label = "Condition on the parentel polygenic risk score?",
           # choices = c("Risk reduction", "Conditional", "Family History"), inline = T
           # choiceValues = c("Risk reduction", "Conditional", "Family History"),
           # choiceNames = c("No conditioning", "Conditional on the parents' polygenic risk score", "Conditional on family disease status"),
           choiceValues = c("Risk reduction", "Conditional"),
           choiceNames = c("No", "Yes"), inline = T,
         ),
         slider_and_numeric(ifelse(num==1, "qf2", "qm2"), 
                            sprintf("Parent %d polygenic risk score percentile:", num), 
                            0.01, 0.99, 0.01, 0.5, paste("For example, if this parameter equal 0.95, the PRS of the father is at the top 5% of the distribution of the PRS in the population.")),
         hr(),
         # Parent 1 Siblings (Aunts/Uncles)
         h5(sprintf("Parent %d Siblings", num)),
         # div(class = "d-flex gap-2",
         layout_column_wrap(width=1/2,
                            # numericInput("sib_p1_sick", "Sick", value = 0, min = 0, width = "50%"),
                            # numericInput("sib_p1_healthy", "Healthy", value = 0, min = 0, width = "50%")
                            sliderInput(sprintf("sib_p%d_n", num),
                                        label = "Number of sibling with known disease status:",
                                        min = 0,
                                        max = 20,
                                        value = 0,
                                        step = 1),
                            sliderInput(sprintf("sib_p%d_sick", num),
                                        label = "Number of affected siblings:",
                                        min = 0,
                                        max = 20,
                                        value = 0,
                                        step = 1)
         )
       ))
}

off <- card(
  card_header("Existing Offspring", class = "bg-secondary text-white"),
  card_body(h6("Siblings of the Embryo"),
            # p(class="text-muted small", "Previous children of this couple"),
            
            # div(class = "d-flex gap-2 justify-content-center",
            layout_column_wrap(width=1/2,
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
                                           step = 1)
                               # numericInput("sib_self_sick", "Sick", value = 0, min = 0, width = "50%"),
                               # numericInput("sib_self_healthy", "Healthy", value = 0, min = 0, width = "50%")
            )
  )
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
  layout_columns(col_widths = c(3, 9), fill=F, 
                 div(
                   embryo_parameters, 
                   disease_parameters),
                   # r2_adj),
                 div(card(card_header("Family information", style = "text-align: center;", class = "bg-primary"),
                          layout_column_wrap(width=1/2,
                                             side(1),
                                             side(2)),
                          off),
                     # layout_column_wrap(width=1/2,
                     # class = "sticky-lg-bottom sticky-cue z-3",
                     card(card_header("Results", class = "bg-success text-white"),  
                          htmlOutput("summary"),
                          class = "sticky-lg-bottom sticky-cue z-3 opacity-75",
                          style = "text-align: center;"))),#),
  disclamir_and_date_text
)

calc_two_traits <- div(class = "well", fluidRow(column(4,
                                                       slider_and_numeric("N_2", "Number of embryos:", 2, 10, 1, 5, "The number of embryos available for selection."),
                                                       slider_and_numeric("rho", '$\\rho$, the genetic correlation between the diseases:', -0.99, 0.99, 0.01, 0, "The genetic correlation between the two diseases."),
                                                       slider_and_numeric("samples_2", "Number of monte carlo draws:", 100000, 500000, 1000, 100000, "The number of simulations. Higher number will give a more accurate estimate, but might take longer to run.")),
                                                column(4, 
                                                       slider_and_numeric("r2_1", "PRS accuracy ($r^2 ~ \\text{disease 1}$):", 0.01, 1, 0.001, 0.05, "The proportion of the variance in the liability of the first disease explained by the polygenic risk score. It is a measure of the accuracy of the score. Typically in the range 0.05-0.1."),
                                                       slider_and_numeric("r2_2", "PRS accuracy ($r^2 ~ \\text{disease 2}$):", 0.01, 1, 0.001, 0.05, "The proportion of the variance in the liability of the second disease explained by the polygenic risk score. It is a measure of the accuracy of the score. Typically in the range 0.05-0.1."),
                                                       fluidRow(column(8, offset = 2, htmlOutput("two_traits"), align = "center"))),
                                                column(4, 
                                                       slider_and_numeric("K_1", "Prevalence of disease 1:", 0.001, 0.3, unique(round(exp(seq(log(0.001), log(0.3), length = 500)), digits = 4)), 0.001, "How prevalent is the first disease in the population? 0.01 means that 1% of the population have the disease, and 0.2 means that 20% of the population have the disease."),
                                                       slider_and_numeric("K_2", "Prevalence of disease 2:", 0.001, 0.3, unique(round(exp(seq(log(0.001), log(0.3), length = 500)), digits = 4)), 0.001, "How prevalent is the second disease in the population? 0.01 means that 1% of the population have the disease, and 0.2 means that 20% of the population have the disease."))),
                       disclamir_and_date_text)

addResourcePath(prefix = "Images", directoryPath = "Images/")

ui <- page_navbar(
  header=withMathJax(tags$head(
  HTML("<script type='text/x-mathjax-config' >
    MathJax.Hub.Config({
    tex2jax: {inlineMath: [['$','$']]}
    });
    </script>"),
  tags$script(HTML(".sticky-cue {
      border-top: 2px solid #dee2e6;
      box-shadow: 0 -4px 8px;
      transition: transform 0.2s ease, padding 0.2s ease;
    }
      .highlight { background:#fff9c4 !important; transition:background 1.5s; }
    .msg { color:#d32f2f; font-size:0.9rem; margin-top:4px;}")),
  tags$script(HTML("Shiny.addCustomMessageHandler('removeClassAfter', function(msg) {
        var el = document.getElementById(msg.id);
        if (!el) return;
        setTimeout(function(){ el.classList.remove(msg.cls); }, msg.delay);
      });"))), 
  shinyjs::useShinyjs()),
  title=tags$a("", uiOutput("logo_ui")),
  nav_panel("About", about_panel),
  nav_panel("Calculator", calc_panel),
  nav_panel("Plot", plot_panel),
  nav_spacer(),
  nav_item(input_dark_mode(id = "dark_mode")),
)
