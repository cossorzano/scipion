# coding: latin-1
# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano
# *
# * Kinestat Pharma
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'coss@cnb.csic.es'
# *
# **************************************************************************
"""
Bibtex string file for Simple package.
"""

_bibtexStr = """

@Article{CHMPEWP56095,
  Title                    = {Guideline on the investigation of drug interactions},
  Author                   = {European Medicines Agency Committee for Human Medicinal Products},
  Journal                  = {CHMP/EWP/560/95},
  Year                     = {2012},
  Doi                      = {http://www.ema.europa.eu/docs/en_GB/document_library/Scientific_guideline/2012/07/WC500129606.pdf},
  Url                      = {http://www.ema.europa.eu/docs/en_GB/document_library/Scientific_guideline/2012/07/WC500129606.pdf}
}

@Article{Fahmi2009,
  Title                    = {Comparison of Different Algorithms for Predicting Clinical Drug-Drug Interactions, Based on the Use of CYP3A4 in Vitro Data: Predictions of Compounds as Precipitants of Interaction},
  Author                   = {Fahmi, O.A. et al},
  Journal                  = {Drug metabolism and disposition},
  Pages                    = {1658-1666},
  Volume                   = {37},
  Year                     = {2009},
  Doi                      = {http://dx.doi.org/10.1124/dmd.108.026252},
  Url                      = {http://dmd.aspetjournals.org/content/37/8/1658}
}

@Article{Kanamitsu2000,
  Title                    = {Quantitative prediction of in vivo drug-drug interactions from in vitro data based on physiological pharmacokinetics: use of maximum unbound concentration of inhibitor at the inlet to the liver},
  Author                   = {Kanamitsu, S. and Ito, K. and Sugiyama, Y.},
  Journal                  = {Pharmaceutical research},
  Pages                    = {336-343},
  Volume                   = {17},
  Year                     = {2000},
  Doi                      = {http://dx.doi.org/10.1023/A:1007572803027},
  Url                      = {http://link.springer.com/article/10.1023/A:1007572803027}
}

@Article{Rostami2004,
  Title                    = {In silico simulations to assess the in vivo consequences of in vitro metabolic drug?drug interactions},
  Author                   = {Rostami-Hodjegan, A. and Tucker, G. T.},
  Journal                  = {Drug Discovery Today},
  Year                     = {2004},
  Pages                    = {441-448},
  Volume                   = {4},
  Doi                      = {http://dx.doi.org/10.1016/j.ddtec.2004.10.002},
  Url                      = {http://www.sciencedirect.com/science/article/pii/S174067490400037X}
}

@Article{Spiess2010,
  Title                    = {An evaluation of R2 as an inadequate measure for nonlinear models in pharmacological and biochemical research: a Monte Carlo approach.},
  Author                   = {Spiess, A. and Neumeyer, N.},
  Journal                  = {BMC Pharmacology},
  Year                     = {2010},
  Pages                    = {6},
  Volume                   = {10},
  Doi                      = {http://dx.doi.org/10.1186/1471-2210-10-6},
  Url                      = {http://dx.doi.org/10.1186/1471-2210-10-6}
}

@Article{Yang2007a,
  Title                    = {Misuse of the Well-Stirred Model of Hepatic Drug Clearance.},
  Author                   = {Yang, J. and Jamei, M. and Yeo, K. R. and Rostami-Hodjegan, A. and Tucker, G. T.},
  Journal                  = {Drug Metabolism and disposition},
  Year                     = {2007},
  Pages                    = {501-502},
  Volume                   = {35},
  Doi                      = {http://dx.doi.org/10.1124/dmd.106.013359},
  Url                      = {http://dx.doi.org/10.1124/dmd.106.013359}
}

@Article{Yang2007b,
  Title                    = {Prediction of Intestinal First-Pass Drug Metabolism.},
  Author                   = {Yang, J. and Jamei, M. and Yeo, K. R. and Tucker, G. T. and Rostami-Hodjegan, A.},
  Journal                  = {Current Drug Metabolism},
  Year                     = {2007},
  Pages                    = {676-684},
  Volume                   = {8},
  Doi                      = {http://dx.doi.org/10.2174/138920007782109733},
  Url                      = {https://www.ncbi.nlm.nih.gov/labs/articles/17979655/}
}

"""

from pyworkflow.utils import parseBibTex

_bibtex = parseBibTex(_bibtexStr)
