20 Sep 2010

Changes since 22 Sep 2008 version:

Test version (printed output is not beautiful) of
madcal.m to perform automatic normalization based on results from IR-MAD
(irmadenvi, irmadss or regirmad) added.

freadenvit.m replaces freadhdr.m which in turn is included in its original
form. freadenvit transposes individual images.

pool.m and poolw.m now do weighted and not simple pooling.

Small updates.

------------------------------------------------------------------

This zip file contains Matlab programs to perform

    - principal component analysis, PCA
    - (weighted) maximum autocorrelation factor analysis, (W)MAF 
    - canonical correlation analysis, CCA
    - (regularized, iteratively reweighted) multivariate alteration detection,
      (IR-)MAD
    - IR-MAD in scale space
    - automatic normalization
    
Since these programs have all the data in memory the analysis should
actually be carried out by means of svd/gsvd; this would improve
numerical stability. However, eig and eigs are used.

Many programs require a small ENVI inspired header file which looks like
this (write the header in any editor)

samples = 150
lines   = 100
bands   = 4
data type = 4
interleave = bsq

where the number of blanks before the equal signs are essential.
Also, no blank after 'bsq' is essential.  Full ENVI headers can be read
but not written.  Headers like these can be written by fwritehdr.

'data type' means
1           uint8
2           int16
3           int32
4           float32
5           float64
6           uint16
7           uint32
8           int64
9           uint64

Presently only 'interleave = bsq' is allowed (no blanks after "bsq"!).

Similar software (with more functionality) written in IDL/ENVI for any-size images is
available from Mort Canty's homepage, http://mcanty.homepage.t-online.de/,
see also M. J. Canty (2010).


Enjoy,

Allan Aasbjerg Nielsen, Ph.D., M.Sc.

Technical University of Denmark
DTU Space - National Space Institute
presently located at
DTU Informatics - Informatics and Mathematical Modelling
Richard Petersens Plads, Building 321
DK-2800 Lyngby, Denmark
aa@space.dtu.dk
http://www.imm.dtu.dk/~aa


References
----------

A. A. Green, M. Berman, P. Switzer, and M. D. Craig (1988).
A Transformation for Ordering Multispectral Data in Terms of Image
Quality with Implications for Noise Removal.
IEEE Transactions on Geoscience and Remote Sensing, 26(1):65-74.

A. A. Nielsen, K. Conradsen and J. J. Simpson (1998).
Multivariate Alteration Detection (MAD) and MAF Postprocessing in
Multispectral, Bitemporal Image Data: New Approaches to Change Detection
Studies.
Remote Sensing of Environment, 64(1):1-19, Elsevier.

T. W. Anderson (2003).
An Introduction to Multivariate Statistics. Third Edition. Wiley.

M. J. Canty, A. A. Nielsen and M. Schmidt (2004).
Automatic Radiometric Normalization of Multitemporal Satellite Imagery.
Remote Sensing of Environment, 91(3-4):441-451.

M. J. Canty and A. A. Nielsen (2004).
Unsupervised Classification of Changes in Multispectral Satellite Imagery.
Proceedings of SPIE, Image and Signal Processing for Remote Sensing X,
5573:356-363.
Maspalomas, Spain.

A. A. Nielsen (2005).
Regularisation in Multi- and Hyperspectral Remote Sensing Change Detection.
6th Geomatic Week Conference, CD-ROM Proceedings,
Barcelona, Spain.

A. A. Nielsen (2005).
An Iterative Extension to the MAD Transformation for Change Detection
in Multi- and Hyper-spectral Remote Sensing Data.
4th EARSeL Workshop on Imaging Spectroscopy, pp. 485-495,
Warsaw, Poland, 27-29 April 2005.

A. A. Nielsen and M. J. Canty (2005).
Multi- and Hyperspectral Remote Sensing Change Detection with Generalized
Difference Images by the IR-MAD Method.
MultiTemp2005, Biloxi, Mississippi, USA.

A. A. Nielsen (2007).
The Regularized Iteratively Reweighted MAD Method for Change Detection
in Multi- and Hyperspectral Data.
IEEE Transactions on Image Processing, 16(2):463-478.

M. J. Canty and A. A. Nielsen (2008).
Automatic Radiometric Normalization of Multitemporal Satellite Imagery
with the Iteratively Re-weighted MAD Transformation.
Remote Sensing of Environment, 112(3):1025-1036.

M. J. Canty (2010).
Image Analysis, Classification and Change Detection in Remote Sensing,
With Algorithms for ENVI/IDL. 2nd revised edition.
Taylor & Francis, CRC Press.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Disclaimer:
%
% No guarantees of performance accompany this software,
% nor is any responsibility assumed on the part of the contributor(s).
% The software has been tested extensively and every effort has been
% made to insure its reliability.
%
% This software is provided by DTU and the contributor(s) 'as is' and
% any express or implied warranties, including, but not limited to, the
% implied warranties of merchantability and fitness for a particular purpose
% are disclaimed.  In no event shall DTU/IMM or the contributor(s) be liable
% for any direct, indirect, incidental, special, exemplary, or consequential
% damages (including, but not limited to, procurement of substitute goods
% or services; loss of use, data, or profits; or business interruption)
% however caused and on any theory of liability, whether in contract, strict
% liability, or tort (including negligence or otherwise) arising in any way
% out of the use of this software, even if advised of the possibility of
% such damage.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
