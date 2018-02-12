""" Create Tex Class


1) Create a Tex File for Calibration

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""        
        
class CreateTex(object):
    @classmethod
    def writeHeader(cls, f, title):
        """Write the header to webpage file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        
        s = '''%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     _               ____  
%    | |        /\   |  _ \ 
%    | |       /  \  | |_) |
%    | |      / /\ \ |  _ < 
%    | |____ / ____ \| |_) |
%    |______/_/    \_\____/ 
%
% _______ ______ __  __ _____  _            _______ ______ 
%|__   __|  ____|  \/  |  __ \| |        /\|__   __|  ____|
%   | |  | |__  | \  / | |__) | |       /  \  | |  | |__   
%   | |  |  __| | |\/| |  ___/| |      / /\ \ | |  |  __|  
%   | |  | |____| |  | | |    | |____ / ____ \| |  | |____ 
%   |_|  |______|_|  |_|_|    |______/_/    \_\_|  |______|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DONT CHANGE ANYTHING BEFORE THE "TITLE" SECTION.%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\\documentclass{article} % Especially this!

\\usepackage[english]{babel}
\\usepackage[utf8]{inputenc}
\\usepackage[margin=1.5in]{geometry}
\\usepackage{amsmath}
\\usepackage{amsthm}
\\usepackage{amsfonts}
\\usepackage{amssymb}
\\usepackage[usenames,dvipsnames]{xcolor}
\\usepackage{graphicx}
\\usepackage{subfig}
\\usepackage[siunitx]{circuitikz}
\\usepackage{tikz}
\\usepackage[colorinlistoftodos, color=orange!50]{todonotes}
\\usepackage{hyperref}
\\usepackage[numbers, square]{natbib}
\\usepackage{fancybox}
\\usepackage{epsfig}
\\usepackage{soul}
\\usepackage[framemethod=tikz]{mdframed}
\\usepackage[shortlabels]{enumitem}
\\usepackage[version=4]{mhchem}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____ _    _  _____ _______ ____  __  __ 
%  / ____| |  | |/ ____|__   __/ __ \|  \/  |
% | |    | |  | | (___    | | | |  | | \  / |
% | |    | |  | |\___ \   | | | |  | | |\/| |
% | |____| |__| |____) |  | | | |__| | |  | |
%  \_____|\____/|_____/   |_|  \____/|_|  |_|
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  _____ ____  __  __ __  __          _   _ _____   _____ 
% / ____/ __ \|  \/  |  \/  |   /\   | \ | |  __ \ / ____|
%| |   | |  | | \  / | \  / |  /  \  |  \| | |  | | (___  
%| |   | |  | | |\/| | |\/| | / /\ \ | . ` | |  | |\___ \ 
%| |___| |__| | |  | | |  | |/ ____ \| |\  | |__| |____) |
% \_____\____/|_|  |_|_|  |_/_/    \_\_| \_|_____/|_____/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SYNTAX FOR NEW COMMANDS:
%\\newcommand{\\new}{Old command or text}

% EXAMPLE:

\\newcommand{\\blah}{blah blah blah \\dots}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  _______ ______          _____ _    _ ______ _____  	%
% |__   __|  ____|   /\   / ____| |  | |  ____|  __ \ 	%
%    | |  | |__     /  \ | |    | |__| | |__  | |__) |	%
%    | |  |  __|   / /\ \| |    |  __  |  __| |  _  / 	%
%    | |  | |____ / ____ \ |____| |  | | |____| | \ \ 	%
%    |_|  |______/_/    \_\_____|_|  |_|______|_|  \_\	%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%														%
% 			COMMANDS				SUMMARY				%
% \\clarity{points}{comment} >>> "Clarity of Writing"	%
% \\other{points}{comment}	>>> "Other"					%
% \\spelling{comment}		>>> "Spelling"				%
% \\units{comment}			>>> "Units"					%
% \\english{comment}			>>> "Syntax and Grammer"	%
% \\source{comment}			>>> "Sources"				%
% \\concept{comment}			>>> "Concept"				%
% \\missing{comment}			>>> "Missing Content"		%
% \\maths{comment}			>>> "Math"					%
% \\terms{comment}			>>> "Science Terms"			%
%														%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\\setlength{\\marginparwidth}{3.4cm}


% NEW COUNTERS
\\newcounter{points}
\\setcounter{points}{100}
\\newcounter{spelling}
\\newcounter{english}
\\newcounter{units}
\\newcounter{other}
\\newcounter{source}
\\newcounter{concept}
\\newcounter{missing}
\\newcounter{math}
\\newcounter{terms}
\\newcounter{clarity}

% COMMANDS

\\definecolor{myblue}{rgb}{0.668, 0.805, 0.929}
\\newcommand{\\hlb}[2][myblue]{ {\\sethlcolor{#1} \\hl{#2}} }

\\newcommand{\\clarity}[2]{\\todo[color=CornflowerBlue!50]{CLARITY of WRITING(#1) #2}\\addtocounter{points}{#1}
\\addtocounter{clarity}{#1}}

\\newcommand{\\other}[2]{\\todo{OTHER(#1) #2} \\addtocounter{points}{#1} \\addtocounter{other}{#1}}

\\newcommand{\\spelling}{\\todo[color=CornflowerBlue!50]{SPELLING (-1)} \\addtocounter{points}{-1}
\\addtocounter{spelling}{-1}}
\\newcommand{\\units}{\\todo{UNITS (-1)} \\addtocounter{points}{-1}
\\addtocounter{units}{-1}}

\\newcommand{\\english}{\\todo[color=CornflowerBlue!50]{SYNTAX and GRAMMAR (-1)} \\addtocounter{points}{-1}
\\addtocounter{english}{-1}}

\\newcommand{\\source}{\\todo{SOURCE(S) (-2)} \\addtocounter{points}{-2}
\\addtocounter{source}{-2}}
\\newcommand{\\concept}{\\todo{CONCEPT (-2)} \\addtocounter{points}{-2}
\\addtocounter{concept}{-2}}

\\newcommand{\\missing}[2]{\\todo{MISSING CONTENT (#1) #2} \\addtocounter{points}{#1}
\\addtocounter{missing}{#1}}

\\newcommand{\\maths}{\\todo{MATH (-1)} \\addtocounter{points}{-1}
\\addtocounter{math}{-1}}
\\newcommand{\\terms}{\\todo[color=CornflowerBlue!50]{SCIENCE TERMS (-1)} \\addtocounter{points}{-1}
\\addtocounter{terms}{-1}}


\\newcommand{\\summary}[1]{
\\begin{mdframed}[nobreak=true]
\\begin{minipage}{\\textwidth}
\\vspace{0.5cm}
\\begin{center}
\\Large{Grade Summary} \\hrule 
\\end{center} \\vspace{0.5cm}
General Comments: #1

\\vspace{0.5cm}
Possible Points \\dotfill 100 \\\\
Points Lost (Science Terms) \\dotfill \\theterms \\\\
Points Lost (Syntax and Grammar) \\dotfill \\theenglish \\\\
Points Lost (Spelling) \\dotfill \\thespelling \\\\
Points Lost (Units) \\dotfill \\theunits \\\\
Points Lost (Math) \\dotfill \\themath \\\\
Points Lost (Sources) \\dotfill \\thesource \\\\
Points Lost (Concept) \\dotfill \\theconcept \\\\
Points Lost (Missing Content) \\dotfill \\themissing \\\\
Points Lost (Clarity of Writing) \\dotfill \\theclarity \\\\
Other \\dotfill \\theother \\\\[0.5cm]
\\begin{center}
\\large{\\textbf{Grade:} \\fbox{\\thepoints}}
\\end{center}
\\end{minipage}
\\end{mdframed}}

%#########################################################

%To use symbols for footnotes
\\renewcommand*{\\thefootnote}{\\fnsymbol{footnote}}
%To change footnotes back to numbers uncomment the following line
%\\renewcommand*{\\thefootnote}{\\arabic{footnote}}

% Enable this command to adjust line spacing for inline math equations.
% \\everymath{\\displaystyle}

% _______ _____ _______ _      ______ 
%|__   __|_   _|__   __| |    |  ____|
%   | |    | |    | |  | |    | |__   
%   | |    | |    | |  | |    |  __|  
%   | |   _| |_   | |  | |____| |____ 
%   |_|  |_____|  |_|  |______|______|
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\\title{
\\rule{\\linewidth}{0.5pt} \\\\[6pt] 
\\LARGE VIRUS Spectrograph Laboratory Calibration Report \\\\
\\Large FILL IN CAMERA ID
\\rule{\\linewidth}{2pt}  \\\\[10pt]
}
\\author{Fill in your name}
\\date{\\normalsize Date Reviewed}

\\begin{document}

\\maketitle
\\noindent


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% _               ____  
%| |        /\   |  _ \ 
%| |       /  \  | |_) |
%| |      / /\ \ |  _ < 
%| |____ / ____ \| |_) |
%|______/_/    \_\____/ 
%%%%%%%%%%%%%%%%%%%%%%%%
%  _____ _______       _____ _______ _____ 
% / ____|__   __|/\   |  __ \__   __/ ____|
%| (___    | |  /  \  | |__) | | | | (___  
% \___ \   | | / /\ \ |  _  /  | |  \___ \ 
% ____) |  | |/ ____ \| | \ \  | |  ____) |
%|_____/   |_/_/    \_\_|  \_\ |_| |_____/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% _    _ ______ _____  ______ 
%| |  | |  ____|  __ \|  ____|
%| |__| | |__  | |__) | |__   
%|  __  |  __| |  _  /|  __|  
%| |  | | |____| | \ \| |____ 
%|_|  |_|______|_|  \_\______|
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

'''

        f.write(s)
    
        
    @classmethod
    def writeObsSummary(cls, f, A, B):
        """Write rows to webpage file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        s = '''
\\section{Overall Summary}
 
% MAKE COMMENTS HERE'''
        t = '''
        
\\begin{table}[H]
\\centering
\\caption{Measured Bias, Gain, and Readnoise}
\\label{basictable}
\\begin{tabular}{|l|l|l|l|}
\\hline
 Amplifier &  Overscan &  Gain (e-/ADU) &  Readnoise (e-) \\\\ \\hline
 %s &  %0.2f & %0.2f & %0.2f   \\\\ \\hline
 %s &  %0.2f & %0.2f & %0.2f   \\\\ \\hline
 %s &  %0.2f & %0.2f & %0.2f   \\\\ \\hline
 %s &  %0.2f & %0.2f & %0.2f   \\\\ \\hline
\\end{tabular}
\\end{table}''' %(tuple(A))
        u = '''
\\begin{table}[H]
\\centering
\\caption{Dark Levels}
\\label{basictable}
\\begin{tabular}{|l|l|l|l|}
\\hline
 Amplifier &  Dark Current & Dark Current &  Dark Current \\\\ 
  &   (DN/pix/sec) &  (e-/pix/sec) &  (e-/pix/600sec) \\\\ \\hline
 %s &  %0.2f & %0.2f & %0.2f   \\\\ \\hline
 %s &  %0.2f & %0.2f & %0.2f   \\\\ \\hline
 %s &  %0.2f & %0.2f & %0.2f   \\\\ \\hline
 %s &  %0.2f & %0.2f & %0.2f   \\\\ \\hline
\\end{tabular}
\\end{table}

''' %(tuple(B))
        f.write(s+t+u)
        f.flush()

    @classmethod
    def writeFigure(cls, f, A):
        """Write rows to webpage file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        u = '''

\\begin{figure}[H]
\\makebox[\\textwidth][c]{\\includegraphics[width=1.5\\textwidth]{%s}}
\\label{testtable}
\\caption{%s}
\\end{figure}


''' % (tuple(A))
        f.write(u)
        f.flush()

    @classmethod
    def writeImageSummary(cls, f, A):
        """Write rows to webpage file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        s = '''
\\section{%s}
''' % A

        t = '''
% MAKE COMMENTS HERE'''
        f.write(s+t)
        f.flush()


    @classmethod 
    def writeEnding(cls, f):
        """Write ending to webpage file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        s = []
        s.append('\\end{document}')
        f.write('\n'.join(s) + "\n")
        f.flush()
        

    