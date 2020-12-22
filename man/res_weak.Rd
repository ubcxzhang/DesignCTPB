\name{res_weak}
\alias{res_weak}
\alias{plot_alpha}
\alias{plot_power}
\alias{opt_alpha}
\alias{opt_power}
\alias{opt_r_split}
\docType{data}
\title{
The clinical trial design for weak biomarker effect
}
\description{
The dataset contains the 3-D rotable figure for the simulated nested two-subset design(weak biomarker effect), from which the optimal choice of each population's proportion and the optimal alpha split can be abtianed. 
}
\usage{data("res_weak")}
\format{
  The format is:
List of 5
 $ plot_power     :List of 8
  ..$ x            :List of 7
  .. ..$ visdat     :List of 1
  .. .. ..$ 7a517404d85b:function ()  
  .. ..$ cur_data   : chr "7a517404d85b"
  .. ..$ attrs      :List of 2
  .. .. ..$ 7a517404d85b:List of 6
  .. .. .. ..$ x           : num [1:101] 0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 ...
  .. .. .. ..$ y           : num [1:101] 0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 ...
  .. .. .. ..$ z           : num [1:101, 1:101] NA NA NA NA NA NA NA NA NA NA ...
  .. .. .. ..$ alpha_stroke: num 1
  .. .. .. ..$ sizes       : num [1:2] 10 100
  .. .. .. ..$ spans       : num [1:2] 1 20
  .. .. ..$ 7a517404d85b:List of 8
  .. .. .. ..$ x           : num [1:101] 0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 ...
  .. .. .. ..$ y           : num [1:101] 0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 ...
  .. .. .. ..$ z           : num [1:101, 1:101] NA NA NA NA NA NA NA NA NA NA ...
  .. .. .. ..$ alpha_stroke: num 1
  .. .. .. ..$ sizes       : num [1:2] 10 100
  .. .. .. ..$ spans       : num [1:2] 1 20
  .. .. .. ..$ type        : chr "surface"
  .. .. .. ..$ inherit     : logi TRUE
  .. .. .. ..- attr(*, "class")= chr [1:2] "plotly_surface" "list"
  .. ..$ layout     :List of 3
  .. .. ..$ width : NULL
  .. .. ..$ height: NULL
  .. .. ..$ margin:List of 4
  .. .. .. ..$ b: num 40
  .. .. .. ..$ l: num 60
  .. .. .. ..$ t: num 25
  .. .. .. ..$ r: num 10
  .. ..$ source     : chr "A"
  .. ..$ config     :List of 1
  .. .. ..$ showSendToCloud: logi FALSE
  .. ..$ layoutAttrs:List of 1
  .. .. ..$ 7a517404d85b:List of 1
  .. .. .. ..$ scene:List of 4
  .. .. .. .. ..$ camera:List of 1
  .. .. .. .. .. ..$ eye:List of 3
  .. .. .. .. .. .. ..$ x: num 2
  .. .. .. .. .. .. ..$ y: num -1
  .. .. .. .. .. .. ..$ z: num 0.34
  .. .. .. .. ..$ xaxis :List of 1
  .. .. .. .. .. ..$ title: chr "r2"
  .. .. .. .. ..$ yaxis :List of 1
  .. .. .. .. .. ..$ title: chr "r3"
  .. .. .. .. ..$ zaxis :List of 1
  .. .. .. .. .. ..$ title: chr "Optimal Power "
  .. ..- attr(*, "TOJSON_FUNC")=function (x, ...)  
  ..$ width        : NULL
  ..$ height       : NULL
  ..$ sizingPolicy :List of 6
  .. ..$ defaultWidth : chr "100%"
  .. ..$ defaultHeight: num 400
  .. ..$ padding      : NULL
  .. ..$ viewer       :List of 6
  .. .. ..$ defaultWidth : NULL
  .. .. ..$ defaultHeight: NULL
  .. .. ..$ padding      : NULL
  .. .. ..$ fill         : logi TRUE
  .. .. ..$ suppress     : logi FALSE
  .. .. ..$ paneHeight   : NULL
  .. ..$ browser      :List of 5
  .. .. ..$ defaultWidth : NULL
  .. .. ..$ defaultHeight: NULL
  .. .. ..$ padding      : NULL
  .. .. ..$ fill         : logi TRUE
  .. .. ..$ external     : logi FALSE
  .. ..$ knitr        :List of 3
  .. .. ..$ defaultWidth : NULL
  .. .. ..$ defaultHeight: NULL
  .. .. ..$ figure       : logi TRUE
  ..$ dependencies :List of 5
  .. ..$ :List of 10
  .. .. ..$ name      : chr "typedarray"
  .. .. ..$ version   : chr "0.1"
  .. .. ..$ src       :List of 1
  .. .. .. ..$ file: chr "htmlwidgets/lib/typedarray"
  .. .. ..$ meta      : NULL
  .. .. ..$ script    : chr "typedarray.min.js"
  .. .. ..$ stylesheet: NULL
  .. .. ..$ head      : NULL
  .. .. ..$ attachment: NULL
  .. .. ..$ package   : chr "plotly"
  .. .. ..$ all_files : logi FALSE
  .. .. ..- attr(*, "class")= chr "html_dependency"
  .. ..$ :List of 10
  .. .. ..$ name      : chr "jquery"
  .. .. ..$ version   : chr "1.11.3"
  .. .. ..$ src       :List of 1
  .. .. .. ..$ file: chr "lib/jquery"
  .. .. ..$ meta      : NULL
  .. .. ..$ script    : chr "jquery.min.js"
  .. .. ..$ stylesheet: NULL
  .. .. ..$ head      : NULL
  .. .. ..$ attachment: NULL
  .. .. ..$ package   : chr "crosstalk"
  .. .. ..$ all_files : logi TRUE
  .. .. ..- attr(*, "class")= chr "html_dependency"
  .. ..$ :List of 10
  .. .. ..$ name      : chr "crosstalk"
  .. .. ..$ version   : chr "1.1.0.1"
  .. .. ..$ src       :List of 1
  .. .. .. ..$ file: chr "www"
  .. .. ..$ meta      : NULL
  .. .. ..$ script    : chr "js/crosstalk.min.js"
  .. .. ..$ stylesheet: chr "css/crosstalk.css"
  .. .. ..$ head      : NULL
  .. .. ..$ attachment: NULL
  .. .. ..$ package   : chr "crosstalk"
  .. .. ..$ all_files : logi TRUE
  .. .. ..- attr(*, "class")= chr "html_dependency"
  .. ..$ :List of 10
  .. .. ..$ name      : chr "plotly-htmlwidgets-css"
  .. .. ..$ version   : chr "1.52.2"
  .. .. ..$ src       :List of 1
  .. .. .. ..$ file: chr "htmlwidgets/lib/plotlyjs"
  .. .. ..$ meta      : NULL
  .. .. ..$ script    : NULL
  .. .. ..$ stylesheet: chr "plotly-htmlwidgets.css"
  .. .. ..$ head      : NULL
  .. .. ..$ attachment: NULL
  .. .. ..$ package   : chr "plotly"
  .. .. ..$ all_files : logi FALSE
  .. .. ..- attr(*, "class")= chr "html_dependency"
  .. ..$ :List of 10
  .. .. ..$ name      : chr "plotly-main"
  .. .. ..$ version   : chr "1.52.2"
  .. .. ..$ src       :List of 1
  .. .. .. ..$ file: chr "htmlwidgets/lib/plotlyjs"
  .. .. ..$ meta      : NULL
  .. .. ..$ script    : chr "plotly-latest.min.js"
  .. .. ..$ stylesheet: NULL
  .. .. ..$ head      : NULL
  .. .. ..$ attachment: NULL
  .. .. ..$ package   : chr "plotly"
  .. .. ..$ all_files : logi FALSE
  .. .. ..- attr(*, "class")= chr "html_dependency"
  ..$ elementId    : NULL
  ..$ preRenderHook:function (p, registerFrames = TRUE)  
  ..$ jsHooks      : list()
  ..- attr(*, "class")= chr [1:2] "plotly" "htmlwidget"
  ..- attr(*, "package")= chr "plotly"
 $ plot_alpha     :List of 8
  ..$ x            :List of 7
  .. ..$ visdat     :List of 4
  .. .. ..$ 7a513a9d3956:function ()  
  .. .. ..$ 7a51728f5bd9:function ()  
  .. .. ..$ 7a51641c67b4:function ()  
  .. .. ..$ 7a51307a50d6:function ()  
  .. ..$ cur_data   : chr "7a51307a50d6"
  .. ..$ attrs      :List of 7
  .. .. ..$ 7a513a9d3956:List of 3
  .. .. .. ..$ alpha_stroke: num 1
  .. .. .. ..$ sizes       : num [1:2] 10 100
  .. .. .. ..$ spans       : num [1:2] 1 20
  .. .. ..$ 7a513a9d3956:List of 8
  .. .. .. ..$ alpha_stroke: num 1
  .. .. .. ..$ sizes       : num [1:2] 10 100
  .. .. .. ..$ spans       : num [1:2] 1 20
  .. .. .. ..$ z           : num [1:101, 1:101] NA NA NA NA NA NA NA NA NA NA ...
  .. .. .. ..$ type        : chr "surface"
  .. .. .. ..$ x           : num [1:101] 0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 ...
  .. .. .. ..$ y           : num [1:101] 0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 ...
  .. .. .. ..$ inherit     : logi TRUE
  .. .. .. ..- attr(*, "class")= chr [1:2] "plotly_surface" "list"
  .. .. ..$ 7a51728f5bd9:List of 12
  .. .. .. ..$ alpha_stroke: num 1
  .. .. .. ..$ sizes       : num [1:2] 10 100
  .. .. .. ..$ spans       : num [1:2] 1 20
  .. .. .. ..$ x           :Class 'formula'  language ~r2
  .. .. .. .. .. ..- attr(*, ".Environment")=<environment: 0x5597a2e01ef0> 
  .. .. .. ..$ y           :Class 'formula'  language ~r3
  .. .. .. .. .. ..- attr(*, ".Environment")=<environment: 0x5597a2e01ef0> 
  .. .. .. ..$ z           :Class 'formula'  language ~alpha1
  .. .. .. .. .. ..- attr(*, ".Environment")=<environment: 0x5597a2e01ef0> 
  .. .. .. ..$ type        : chr "scatter3d"
  .. .. .. ..$ mode        : chr "markers"
  .. .. .. ..$ size        : num 2
  .. .. .. ..$ symbol      : num 0
  .. .. .. ..$ name        : chr "alpha1"
  .. .. .. ..$ inherit     : logi TRUE
  .. .. ..$ 7a51728f5bd9:List of 9
  .. .. .. ..$ alpha_stroke: num 1
  .. .. .. ..$ sizes       : num [1:2] 10 100
  .. .. .. ..$ spans       : num [1:2] 1 20
  .. .. .. ..$ z           : num [1:101, 1:101] NA NA NA NA NA NA NA NA NA NA ...
  .. .. .. ..$ type        : chr "surface"
  .. .. .. ..$ x           : num [1:101] 0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 ...
  .. .. .. ..$ y           : num [1:101] 0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 ...
  .. .. .. ..$ opacity     : num 0.98
  .. .. .. ..$ inherit     : logi TRUE
  .. .. .. ..- attr(*, "class")= chr [1:2] "plotly_surface" "list"
  .. .. ..$ 7a51641c67b4:List of 12
  .. .. .. ..$ alpha_stroke: num 1
  .. .. .. ..$ sizes       : num [1:2] 10 100
  .. .. .. ..$ spans       : num [1:2] 1 20
  .. .. .. ..$ x           :Class 'formula'  language ~r2
  .. .. .. .. .. ..- attr(*, ".Environment")=<environment: 0x5597a2dfeb30> 
  .. .. .. ..$ y           :Class 'formula'  language ~r3
  .. .. .. .. .. ..- attr(*, ".Environment")=<environment: 0x5597a2dfeb30> 
  .. .. .. ..$ z           :Class 'formula'  language ~alpha2
  .. .. .. .. .. ..- attr(*, ".Environment")=<environment: 0x5597a2dfeb30> 
  .. .. .. ..$ type        : chr "scatter3d"
  .. .. .. ..$ mode        : chr "markers"
  .. .. .. ..$ size        : num 2
  .. .. .. ..$ symbol      : num 100
  .. .. .. ..$ name        : chr "alpha2"
  .. .. .. ..$ inherit     : logi TRUE
  .. .. ..$ 7a51641c67b4:List of 9
  .. .. .. ..$ alpha_stroke: num 1
  .. .. .. ..$ sizes       : num [1:2] 10 100
  .. .. .. ..$ spans       : num [1:2] 1 20
  .. .. .. ..$ z           : num [1:101, 1:101] NA NA NA NA NA NA NA NA NA NA ...
  .. .. .. ..$ type        : chr "surface"
  .. .. .. ..$ x           : num [1:101] 0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 ...
  .. .. .. ..$ y           : num [1:101] 0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 ...
  .. .. .. ..$ opacity     : num 0.98
  .. .. .. ..$ inherit     : logi TRUE
  .. .. .. ..- attr(*, "class")= chr [1:2] "plotly_surface" "list"
  .. .. ..$ 7a51307a50d6:List of 12
  .. .. .. ..$ alpha_stroke: num 1
  .. .. .. ..$ sizes       : num [1:2] 10 100
  .. .. .. ..$ spans       : num [1:2] 1 20
  .. .. .. ..$ x           :Class 'formula'  language ~r2
  .. .. .. .. .. ..- attr(*, ".Environment")=<environment: 0x5597a2df57f0> 
  .. .. .. ..$ y           :Class 'formula'  language ~r3
  .. .. .. .. .. ..- attr(*, ".Environment")=<environment: 0x5597a2df57f0> 
  .. .. .. ..$ z           :Class 'formula'  language ~alpha3
  .. .. .. .. .. ..- attr(*, ".Environment")=<environment: 0x5597a2df57f0> 
  .. .. .. ..$ type        : chr "scatter3d"
  .. .. .. ..$ mode        : chr "markers"
  .. .. .. ..$ size        : num 2
  .. .. .. ..$ symbol      : num 200
  .. .. .. ..$ name        : chr "alpha3"
  .. .. .. ..$ inherit     : logi TRUE
  .. ..$ layout     :List of 3
  .. .. ..$ width : NULL
  .. .. ..$ height: NULL
  .. .. ..$ margin:List of 4
  .. .. .. ..$ b: num 40
  .. .. .. ..$ l: num 60
  .. .. .. ..$ t: num 25
  .. .. .. ..$ r: num 10
  .. ..$ source     : chr "A"
  .. ..$ config     :List of 1
  .. .. ..$ showSendToCloud: logi FALSE
  .. ..$ layoutAttrs:List of 1
  .. .. ..$ 7a51307a50d6:List of 1
  .. .. .. ..$ scene:List of 4
  .. .. .. .. ..$ camera:List of 1
  .. .. .. .. .. ..$ eye:List of 3
  .. .. .. .. .. .. ..$ x: num 2
  .. .. .. .. .. .. ..$ y: num -1
  .. .. .. .. .. .. ..$ z: num 0.34
  .. .. .. .. ..$ xaxis :List of 1
  .. .. .. .. .. ..$ title: chr "r2"
  .. .. .. .. ..$ yaxis :List of 1
  .. .. .. .. .. ..$ title: chr "r3"
  .. .. .. .. ..$ zaxis :List of 1
  .. .. .. .. .. ..$ title: chr "Optimal alpha"
  .. ..- attr(*, "TOJSON_FUNC")=function (x, ...)  
  ..$ width        : NULL
  ..$ height       : NULL
  ..$ sizingPolicy :List of 6
  .. ..$ defaultWidth : chr "100%"
  .. ..$ defaultHeight: num 400
  .. ..$ padding      : NULL
  .. ..$ viewer       :List of 6
  .. .. ..$ defaultWidth : NULL
  .. .. ..$ defaultHeight: NULL
  .. .. ..$ padding      : NULL
  .. .. ..$ fill         : logi TRUE
  .. .. ..$ suppress     : logi FALSE
  .. .. ..$ paneHeight   : NULL
  .. ..$ browser      :List of 5
  .. .. ..$ defaultWidth : NULL
  .. .. ..$ defaultHeight: NULL
  .. .. ..$ padding      : NULL
  .. .. ..$ fill         : logi TRUE
  .. .. ..$ external     : logi FALSE
  .. ..$ knitr        :List of 3
  .. .. ..$ defaultWidth : NULL
  .. .. ..$ defaultHeight: NULL
  .. .. ..$ figure       : logi TRUE
  ..$ dependencies :List of 5
  .. ..$ :List of 10
  .. .. ..$ name      : chr "typedarray"
  .. .. ..$ version   : chr "0.1"
  .. .. ..$ src       :List of 1
  .. .. .. ..$ file: chr "htmlwidgets/lib/typedarray"
  .. .. ..$ meta      : NULL
  .. .. ..$ script    : chr "typedarray.min.js"
  .. .. ..$ stylesheet: NULL
  .. .. ..$ head      : NULL
  .. .. ..$ attachment: NULL
  .. .. ..$ package   : chr "plotly"
  .. .. ..$ all_files : logi FALSE
  .. .. ..- attr(*, "class")= chr "html_dependency"
  .. ..$ :List of 10
  .. .. ..$ name      : chr "jquery"
  .. .. ..$ version   : chr "1.11.3"
  .. .. ..$ src       :List of 1
  .. .. .. ..$ file: chr "lib/jquery"
  .. .. ..$ meta      : NULL
  .. .. ..$ script    : chr "jquery.min.js"
  .. .. ..$ stylesheet: NULL
  .. .. ..$ head      : NULL
  .. .. ..$ attachment: NULL
  .. .. ..$ package   : chr "crosstalk"
  .. .. ..$ all_files : logi TRUE
  .. .. ..- attr(*, "class")= chr "html_dependency"
  .. ..$ :List of 10
  .. .. ..$ name      : chr "crosstalk"
  .. .. ..$ version   : chr "1.1.0.1"
  .. .. ..$ src       :List of 1
  .. .. .. ..$ file: chr "www"
  .. .. ..$ meta      : NULL
  .. .. ..$ script    : chr "js/crosstalk.min.js"
  .. .. ..$ stylesheet: chr "css/crosstalk.css"
  .. .. ..$ head      : NULL
  .. .. ..$ attachment: NULL
  .. .. ..$ package   : chr "crosstalk"
  .. .. ..$ all_files : logi TRUE
  .. .. ..- attr(*, "class")= chr "html_dependency"
  .. ..$ :List of 10
  .. .. ..$ name      : chr "plotly-htmlwidgets-css"
  .. .. ..$ version   : chr "1.52.2"
  .. .. ..$ src       :List of 1
  .. .. .. ..$ file: chr "htmlwidgets/lib/plotlyjs"
  .. .. ..$ meta      : NULL
  .. .. ..$ script    : NULL
  .. .. ..$ stylesheet: chr "plotly-htmlwidgets.css"
  .. .. ..$ head      : NULL
  .. .. ..$ attachment: NULL
  .. .. ..$ package   : chr "plotly"
  .. .. ..$ all_files : logi FALSE
  .. .. ..- attr(*, "class")= chr "html_dependency"
  .. ..$ :List of 10
  .. .. ..$ name      : chr "plotly-main"
  .. .. ..$ version   : chr "1.52.2"
  .. .. ..$ src       :List of 1
  .. .. .. ..$ file: chr "htmlwidgets/lib/plotlyjs"
  .. .. ..$ meta      : NULL
  .. .. ..$ script    : chr "plotly-latest.min.js"
  .. .. ..$ stylesheet: NULL
  .. .. ..$ head      : NULL
  .. .. ..$ attachment: NULL
  .. .. ..$ package   : chr "plotly"
  .. .. ..$ all_files : logi FALSE
  .. .. ..- attr(*, "class")= chr "html_dependency"
  ..$ elementId    : NULL
  ..$ preRenderHook:function (p, registerFrames = TRUE)  
  ..$ jsHooks      : list()
  ..- attr(*, "class")= chr [1:2] "plotly" "htmlwidget"
  ..- attr(*, "package")= chr "plotly"
 $ opt_r_split    : Named num [1:2] 0.303 0
  ..- attr(*, "names")= chr [1:2] "r2" "r3"
 $ opt_power      : num 0.788
 $ opt_alpha_split: Named num [1:3] 0.01254 0.00844 0.00599
  ..- attr(*, "names")= chr [1:3] "alpha1" "alpha2" "alpha3"
}
\examples{
data(res_weak)
}
\keyword{datasets}
