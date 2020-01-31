def initMPL():
    from matplotlib import rc

    font = {
    'family'        : 'CMU Serif Roman',
    'serif'         : 'CMU Serif Roman',
    'sans-serif'    : 'CMU Serif Roman',
    'cursive'       : 'DejaVu Sans',
    'monospace'     : 'CMU Typewriter',
}
    rc('font',**font)
    ## for Palatino and other serif fonts use:
    #rc('font',**{'family':'serif','serif':['Palatino']})
    rc('text', usetex=True)
    rc('text.latex', preamble=r"\usepackage{amsmath}\newcommand{\colvec}[1]{\ensuremath{\begin{pmatrix}#1\end{pmatrix}}}")
    x = "'164"
    y = "'171"
    #rc('text.latex', preamble=r'\let\llangle\@undefined\let\rrangle\@undefined\DeclareMathDelimiter{\llangle}{\mathopen}{MnLargeSymbols}{'+x+'}{MnLargeSymbols}{'+x+'}\DeclareMathDelimiter{\rrangle}{\mathclose}{MnLargeSymbols}{'+y+'}{MnLargeSymbols}{'+y+'}')
