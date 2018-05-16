MathJax.Hub.Config({
  extensions: ["MathMenu.js", "MathZoom.js", "tex2jax.js"],
  config: ["MMLorHTML.js"],
  jax: ["input/TeX", "output/HTML-CSS", "output/NativeMML"],
  tex2jax: { inlineMath: [ [ '$', '$' ] ] },
  TeX: { equationNumbers: { autoNumber: "AMS"},
         Macros: { 
            RR:    '{\\mathbf R}', 
            bold: ['\\mathbf{#1}', 1]
                 }       
       };
});
