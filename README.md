# Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap

This repository contains the data, code, and analyses used in Gandal et al., (2017) "Shared molecular neuropathology across major psychiatric disorders parallels polygenic overlap" [doi: 10.1101/040022](https://doi.org/10.1101/040022)

### Organization
1. `code/` contains the code for normalizing individual datasets and performing cross disroder transcriptome analyess
2. `raw_data/` contains individual-level study data for microarray,  RNAseq (counts), and GWAS (sum-stats)
3. `results/` contains figures and table results from the cross disorder analyses
4. `working_data/` contains intermediary files


### Microarray datasets used in this study

<table class=MsoNormalTable border=0 cellspacing=0 cellpadding=0 width=464
 style='width:464.0pt;border-collapse:collapse;mso-yfti-tbllook:1184;
 mso-padding-alt:0in 5.4pt 0in 5.4pt'>
 <tr style='mso-yfti-irow:0;mso-yfti-firstrow:yes;height:16.0pt'>
  <td width=51 rowspan=2 style='width:50.5pt;border:solid black 1.0pt;
  background:#4E5B6F;padding:0in 1.4pt 0in 1.4pt;height:16.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><b><span
  style='font-size:8.0pt;color:white'>Disease<o:p></o:p></span></b></p>
  </td>
  <td width=109 colspan=2 style='width:109.25pt;border-top:solid black 1.0pt;
  border-left:none;border-bottom:none;border-right:solid black 1.0pt;
  background:#4E5B6F;padding:0in 1.4pt 0in 1.4pt;height:16.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><b><span
  style='font-size:8.0pt;color:white'># Samples<o:p></o:p></span></b></p>
  </td>
  <td width=72 rowspan=2 style='width:72.3pt;border:solid black 1.0pt;
  border-left:none;mso-border-left-alt:solid black 1.0pt;background:#4E5B6F;
  padding:0in 1.4pt 0in 1.4pt;height:16.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><b><span
  style='font-size:8.0pt;color:white'>Brain Region<o:p></o:p></span></b></p>
  </td>
  <td width=103 rowspan=2 style='width:102.75pt;border:solid black 1.0pt;
  border-left:none;mso-border-left-alt:solid black 1.0pt;background:#4E5B6F;
  padding:0in 1.4pt 0in 1.4pt;height:16.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><b><span
  style='font-size:8.0pt;color:white'>Platform<o:p></o:p></span></b></p>
  </td>
  <td width=55 rowspan=2 style='width:55.2pt;border:solid black 1.0pt;
  border-left:none;mso-border-left-alt:solid black 1.0pt;background:#4E5B6F;
  padding:0in 1.4pt 0in 1.4pt;height:16.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><b><span
  style='font-size:8.0pt;color:white'>Study<o:p></o:p></span></b></p>
  </td>
  <td width=74 rowspan=2 style='width:74.0pt;border:solid black 1.0pt;
  border-left:none;mso-border-left-alt:solid black 1.0pt;background:#4E5B6F;
  padding:0in 1.4pt 0in 1.4pt;height:16.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><b><span
  style='font-size:8.0pt;color:white'>Dataset<o:p></o:p></span></b></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:1;height:3.05pt'>
  <td width=55 style='width:55.25pt;border:none;border-bottom:solid black 1.0pt;
  background:#4E5B6F;padding:0in 1.4pt 0in 1.4pt;height:3.05pt'>
  <p class=MsoNormal align=center style='text-align:center'><b><span
  style='font-size:8.0pt;color:white'>Cases<o:p></o:p></span></b></p>
  </td>
  <td width=54 style='width:.75in;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:#4E5B6F;
  padding:0in 1.4pt 0in 1.4pt;height:3.05pt'>
  <p class=MsoNormal align=center style='text-align:center'><b><span
  style='font-size:8.0pt;color:white'>Controls<o:p></o:p></span></b></p>
  </td>

 </tr>
 <tr style='mso-yfti-irow:2;height:9.35pt'>
  <td width=51 rowspan=3 style='width:50.5pt;border:solid black 1.0pt;
  border-top:none;padding:0in 1.4pt 0in 1.4pt;height:9.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><b><span
  style='font-size:8.0pt;color:black'>ASD<o:p></o:p></span></b></p>
  </td>
  <td width=55 style='width:55.25pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;background:
  #D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:9.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>29<o:p></o:p></span></p>
  </td>
  <td width=54 style='width:.75in;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:#D9D9D9;
  padding:0in 1.4pt 0in 1.4pt;height:9.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>29<o:p></o:p></span></p>
  </td>
  <td width=72 style='width:72.3pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:#D9D9D9;
  padding:0in 1.4pt 0in 1.4pt;height:9.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>BA9, BA41<o:p></o:p></span></p>
  </td>
  <td width=103 style='width:102.75pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;background:
  #D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:9.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Illumina Ref8 v3<o:p></o:p></span></p>
  </td>
  <td width=55 style='width:55.2pt;border:none;border-right:solid black 1.0pt;
  background:#D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:9.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Voineagu<o:p></o:p></span></p>
  </td>
  <td width=74 style='width:74.0pt;border:none;border-right:solid black 1.0pt;
  background:#D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:9.35pt'><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28521">GSE28521</a>
  </td>
 </tr>
 <tr style='mso-yfti-irow:3;height:6.65pt'>
  <td width=55 style='width:55.25pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;background:
  white;padding:0in 1.4pt 0in 1.4pt;height:6.65pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>15<o:p></o:p></span></p>
  </td>
  <td width=54 style='width:.75in;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:white;padding:
  0in 5.4pt 0in 5.4pt;height:6.65pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>18<o:p></o:p></span></p>
  </td>
  <td width=72 style='width:72.3pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:white;padding:
  0in 5.4pt 0in 5.4pt;height:6.65pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>BA9/46<o:p></o:p></span></p>
  </td>
  <td width=103 style='width:102.75pt;border:none;border-bottom:solid black 1.0pt;
  background:white;padding:0in 1.4pt 0in 1.4pt;height:6.65pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Illumina Ref8 v3<o:p></o:p></span></p>
  </td>
  <td width=55 style='width:55.2pt;border:solid black 1.0pt;background:white;
  padding:0in 1.4pt 0in 1.4pt;height:6.65pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Chow<o:p></o:p></span></p>
  </td>
  <td width=74 style='width:74.0pt;border:solid black 1.0pt;border-left:none;
  background:white;padding:0in 1.4pt 0in 1.4pt;height:6.65pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28475">GSE28475</a><o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:4;height:17.0pt'>
  <td width=55 style='width:55.25pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;background:
  #D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:17.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>6<o:p></o:p></span></p>
  </td>
  <td width=54 style='width:.75in;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:#D9D9D9;
  padding:0in 1.4pt 0in 1.4pt;height:17.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>6<o:p></o:p></span></p>
  </td>
  <td width=72 style='width:72.3pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:#D9D9D9;
  padding:0in 1.4pt 0in 1.4pt;height:17.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>BA41/42<o:p></o:p></span></p>
  </td>
  <td width=103 style='width:102.75pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;background:
  #D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:17.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Affy HG-U133 plus2<o:p></o:p></span></p>
  </td>
  <td width=129 colspan=2 style='width:129.2pt;border-top:none;border-left:
  none;border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;
  background:#D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:17.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'><a href="https://www.ncbi.nlm.nih.gov/pubmed/18378158">Garbett et al., 2008</a></span>
  </td>
 </tr>
 <tr style='mso-yfti-irow:5;height:29.0pt'>
  <td width=51 rowspan=6 style='width:50.5pt;border:solid black 1.0pt;
  border-top:none;padding:0in 1.4pt 0in 1.4pt;height:29.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><b><span
  style='font-size:8.0pt;color:black'>SCZ<o:p></o:p></span></b></p>
  </td>
  <td width=55 rowspan=2 style='width:55.25pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;mso-border-left-alt:
  solid black 1.0pt;background:white;padding:0in 1.4pt 0in 1.4pt;height:29.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>51<o:p></o:p></span></p>
  </td>
  <td width=54 rowspan=2 style='width:.75in;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;mso-border-left-alt:
  solid black 1.0pt;background:white;padding:0in 1.4pt 0in 1.4pt;height:29.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>50<o:p></o:p></span></p>
  </td>
  <td width=72 rowspan=2 style='width:72.3pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;mso-border-left-alt:
  solid black 1.0pt;background:white;padding:0in 1.4pt 0in 1.4pt;height:29.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Parietal cortex<o:p></o:p></span></p>
  </td>
  <td width=103 rowspan=2 style='width:102.75pt;border-top:none;border-left:
  none;border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;
  mso-border-left-alt:solid black 1.0pt;background:white;padding:0in 1.4pt 0in 1.4pt;
  height:29.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Affy HuGene 1.0 ST<o:p></o:p></span></p>
  </td>
  <td width=55 rowspan=2 style='width:55.2pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;mso-border-left-alt:
  solid black 1.0pt;background:white;padding:0in 1.4pt 0in 1.4pt;height:29.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Chen<o:p></o:p></span></p>
  </td>
  <td width=74 rowspan=2 style='width:74.0pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;mso-border-left-alt:
  solid black 1.0pt;background:white;padding:0in 1.4pt 0in 1.4pt;height:29.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:#211E1E'><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35978">GSE35978</a><o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:6;height:12.05pt'>
   <tr style='mso-yfti-irow:7;height:3.35pt'>
  <td width=55 style='width:55.25pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;background:
  #D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:3.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>15<o:p></o:p></span></p>
  </td>
  <td width=54 style='width:.75in;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:#D9D9D9;
  padding:0in 1.4pt 0in 1.4pt;height:3.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>19<o:p></o:p></span></p>
  </td>
  <td width=72 style='width:72.3pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:#D9D9D9;
  padding:0in 1.4pt 0in 1.4pt;height:3.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>BA46<o:p></o:p></span></p>
  </td>
  <td width=103 style='width:102.75pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;background:
  #D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:3.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Affy HG-U133 plus2 <o:p></o:p></span></p>
  </td>
  <td width=55 style='width:55.2pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:#D9D9D9;
  padding:0in 1.4pt 0in 1.4pt;height:3.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Lanz<o:p></o:p></span></p>
  </td>
  <td width=74 style='width:74.0pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:#D9D9D9;
  padding:0in 1.4pt 0in 1.4pt;height:3.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53987">GSE53987</a><o:p></o:p></span></p>
  </td>
  </tr>
 <tr style='mso-yfti-irow:8;height:3.35pt'>
  <td width=55 style='width:55.25pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;background:
  white;padding:0in 1.4pt 0in 1.4pt;height:3.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>28<o:p></o:p></span></p>
  </td>
  <td width=54 style='width:.75in;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:white;padding:
  0in 5.4pt 0in 5.4pt;height:3.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>23<o:p></o:p></span></p>
  </td>
  <td width=72 style='width:72.3pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:white;padding:
  0in 5.4pt 0in 5.4pt;height:3.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>BA10<o:p></o:p></span></p>
  </td>
  <td width=103 style='width:102.75pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;background:
  white;padding:0in 1.4pt 0in 1.4pt;height:3.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Affy HG-U133 plus2<o:p></o:p></span></p>
  </td>
  <td width=55 style='width:55.2pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:white;padding:
  0in 5.4pt 0in 5.4pt;height:3.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Maycox<o:p></o:p></span></p>
  </td>
  <td width=74 style='width:74.0pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:white;padding:
  0in 5.4pt 0in 5.4pt;height:3.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17612">GSE17612</a><o:p></o:p></span></p>
  </td>
  
 </tr>
 <tr style='mso-yfti-irow:9;height:3.35pt'>
  <td width=55 style='width:55.25pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;background:
  #D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:3.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>35<o:p></o:p></span></p>
  </td>
  <td width=54 style='width:.75in;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:#D9D9D9;
  padding:0in 1.4pt 0in 1.4pt;height:3.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>34<o:p></o:p></span></p>
  </td>
  <td width=72 style='width:72.3pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:#D9D9D9;
  padding:0in 1.4pt 0in 1.4pt;height:3.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>BA46<o:p></o:p></span></p>
  </td>
  <td width=103 style='width:102.75pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;background:
  #D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:3.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Affy HG-U133A<o:p></o:p></span></p>
  </td>
  <td width=55 style='width:55.2pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:#D9D9D9;
  padding:0in 1.4pt 0in 1.4pt;height:3.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Iwamoto<o:p></o:p></span></p>
  </td>
  <td width=74 style='width:74.0pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:#D9D9D9;
  padding:0in 1.4pt 0in 1.4pt;height:3.35pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12649">GSE12649</a><o:p></o:p></span></p>
  </td>
  
 </tr>
 <tr style='mso-yfti-irow:10;height:15.2pt'>
  <td width=55 style='width:55.25pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;background:
  white;padding:0in 1.4pt 0in 1.4pt;height:15.2pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>30<o:p></o:p></span></p>
  </td>
  <td width=54 style='width:.75in;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:white;padding:
  0in 5.4pt 0in 5.4pt;height:15.2pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>29<o:p></o:p></span></p>
  </td>
  <td width=72 style='width:72.3pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:white;padding:
  0in 5.4pt 0in 5.4pt;height:15.2pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>BA46<o:p></o:p></span></p>
  </td>
  <td width=103 style='width:102.75pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;background:
  white;padding:0in 1.4pt 0in 1.4pt;height:15.2pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Affy HG-U133 plus2<o:p></o:p></span></p>
  </td>
  <td width=55 style='width:55.2pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:white;padding:
  0in 5.4pt 0in 5.4pt;height:15.2pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Narayan<o:p></o:p></span></p>
  </td>
  <td width=74 style='width:74.0pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:white;padding:
  0in 5.4pt 0in 5.4pt;height:15.2pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21138">GSE21138</a><o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:11;height:29.0pt'>
  <td width=51 rowspan=4 style='width:50.5pt;border:solid black 1.0pt;
  border-top:none;padding:0in 1.4pt 0in 1.4pt;height:29.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><b><span
  style='font-size:8.0pt;color:black'>BD<o:p></o:p></span></b></p>
  </td>
  <td width=55 rowspan=2 style='width:55.25pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;mso-border-left-alt:
  solid black 1.0pt;background:#D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:29.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>45<o:p></o:p></span></p>
  </td>
  <td width=54 rowspan=2 style='width:.75in;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;mso-border-left-alt:
  solid black 1.0pt;background:#D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:29.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Included above (50)<o:p></o:p></span></p>
  </td>
  <td width=72 rowspan=2 style='width:72.3pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;mso-border-left-alt:
  solid black 1.0pt;background:#D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:29.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Parietal cortex<o:p></o:p></span></p>
  </td>
  <td width=103 rowspan=2 style='width:102.75pt;border-top:none;border-left:
  none;border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;
  mso-border-left-alt:solid black 1.0pt;background:#D9D9D9;padding:0in 1.4pt 0in 1.4pt;
  height:29.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Affy HuGene 1.0 ST<o:p></o:p></span></p>
  </td>
  <td width=55 rowspan=2 style='width:55.2pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;mso-border-left-alt:
  solid black 1.0pt;background:#D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:29.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Chen<o:p></o:p></span></p>
  </td>
  <td width=74 rowspan=2 style='width:74.0pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;mso-border-left-alt:
  solid black 1.0pt;background:#D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:29.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:#211E1E'><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35978">GSE35978</a><o:p></o:p></span></p>
  </td>
   </tr>
 <tr style='mso-yfti-irow:12;height:12.05pt'>
  
 </tr>
 <tr style='mso-yfti-irow:13;height:4.0pt'>
  <td width=55 style='width:55.25pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;background:
  white;padding:0in 1.4pt 0in 1.4pt;height:4.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>17<o:p></o:p></span></p>
  </td>
  <td width=54 style='width:.75in;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:white;padding:
  0in 5.4pt 0in 5.4pt;height:4.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Included above (19)<o:p></o:p></span></p>
  </td>
  <td width=72 style='width:72.3pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:white;padding:
  0in 5.4pt 0in 5.4pt;height:4.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>BA46<o:p></o:p></span></p>
  </td>
  <td width=103 style='width:102.75pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;background:
  white;padding:0in 1.4pt 0in 1.4pt;height:4.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Affy HG-U133 plus2 <o:p></o:p></span></p>
  </td>
  <td width=55 style='width:55.2pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:white;padding:
  0in 5.4pt 0in 5.4pt;height:4.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Lanz<o:p></o:p></span></p>
  </td>
  <td width=74 style='width:74.0pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:white;padding:
  0in 5.4pt 0in 5.4pt;height:4.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53987">GSE53987</a><o:p></o:p></span></p>
  </td>
  
 </tr>
 <tr style='mso-yfti-irow:14;height:15.2pt'>
  <td width=55 style='width:55.25pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;background:
  #D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:15.2pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>32<o:p></o:p></span></p>
  </td>
  <td width=54 style='width:.75in;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:#D9D9D9;
  padding:0in 1.4pt 0in 1.4pt;height:15.2pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Included above (34)<o:p></o:p></span></p>
  </td>
  <td width=72 style='width:72.3pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:#D9D9D9;
  padding:0in 1.4pt 0in 1.4pt;height:15.2pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>BA46<o:p></o:p></span></p>
  </td>
  <td width=103 style='width:102.75pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;background:
  #D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:15.2pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Affy HG-U133A<o:p></o:p></span></p>
  </td>
  <td width=55 style='width:55.2pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:#D9D9D9;
  padding:0in 1.4pt 0in 1.4pt;height:15.2pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Iwamoto<o:p></o:p></span></p>
  </td>
  <td width=74 style='width:74.0pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:#D9D9D9;
  padding:0in 1.4pt 0in 1.4pt;height:15.2pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12649">GSE12649</a><o:p></o:p></span></p>
  </td>
  
 </tr>
 <tr style='mso-yfti-irow:15;height:5.3pt'>
  <td width=51 rowspan=5 style='width:50.5pt;border:solid black 1.0pt;
  border-top:none;padding:0in 1.4pt 0in 1.4pt;height:5.3pt'>
  <p class=MsoNormal align=center style='text-align:center'><b><span
  style='font-size:8.0pt;color:black'>MDD<o:p></o:p></span></b></p>
  </td>
  <td width=55 style='width:55.25pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;background:
  white;padding:0in 1.4pt 0in 1.4pt;height:5.3pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>17<o:p></o:p></span></p>
  </td>
  <td width=54 style='width:.75in;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:white;padding:
  0in 5.4pt 0in 5.4pt;height:5.3pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Included above (19)<o:p></o:p></span></p>
  </td>
  <td width=72 style='width:72.3pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:white;padding:
  0in 5.4pt 0in 5.4pt;height:5.3pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>BA46<o:p></o:p></span></p>
  </td>
  <td width=103 style='width:102.75pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;background:
  white;padding:0in 1.4pt 0in 1.4pt;height:5.3pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Affy HG-U133 plus2<o:p></o:p></span></p>
  </td>
  <td width=55 style='width:55.2pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:white;padding:
  0in 5.4pt 0in 5.4pt;height:5.3pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Lanz<o:p></o:p></span></p>
  </td>
  <td width=74 style='width:74.0pt;border-top:none;border-left:none;border-bottom:
  solid black 1.0pt;border-right:solid black 1.0pt;background:white;padding:
  0in 5.4pt 0in 5.4pt;height:5.3pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53987">GSE53987</a><o:p></o:p></span></p>
  </td>
  
 </tr>
 <tr style='mso-yfti-irow:16;height:16.0pt'>
  <td width=55 rowspan=4 style='width:55.25pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;mso-border-left-alt:
  solid black 1.0pt;background:#D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:16.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>70<o:p></o:p></span></p>
  </td>
  <td width=54 rowspan=4 style='width:.75in;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;mso-border-left-alt:
  solid black 1.0pt;background:#D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:16.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>70<o:p></o:p></span></p>
  </td>
  <td width=72 rowspan=4 style='width:72.3pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;mso-border-left-alt:
  solid black 1.0pt;background:#D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:16.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>BA9, BA25<o:p></o:p></span></p>
  </td>
  <td width=103 rowspan=4 style='width:102.75pt;border-top:none;border-left:
  none;border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;
  mso-border-left-alt:solid black 1.0pt;background:#D9D9D9;padding:0in 1.4pt 0in 1.4pt;
  height:16.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Affy HG-U133 plus2<o:p></o:p></span></p>
  </td>
  <td width=55 rowspan=4 style='width:55.2pt;border-top:none;border-left:none;
  border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;mso-border-left-alt:
  solid black 1.0pt;background:#D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:16.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Sibille<o:p></o:p></span></p>
  </td>
  <td width=74 valign=bottom style='width:74.0pt;border:none;border-right:solid black 1.0pt;
  background:#D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:16.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54567">GSE54567</a><o:p></o:p></span></p>
  </td>
  
 </tr>
 <tr style='mso-yfti-irow:17;height:16.0pt'>
  <td width=74 valign=bottom style='width:74.0pt;border:none;border-right:solid black 1.0pt;
  background:#D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:16.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54568">GSE54568</a><o:p></o:p></span></p>
  </td>
  
 </tr>
 <tr style='mso-yfti-irow:18;height:16.0pt'>
  <td width=74 valign=bottom style='width:74.0pt;border:none;border-right:solid black 1.0pt;
  background:#D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:16.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54571">GSE54571</a><o:p></o:p></span></p>
  </td>
  
 </tr>
 <tr style='mso-yfti-irow:19;height:4.0pt'>
  <td width=74 valign=bottom style='width:74.0pt;border-top:none;border-left:
  none;border-bottom:solid black 1.0pt;border-right:solid black 1.0pt;
  background:#D9D9D9;padding:0in 1.4pt 0in 1.4pt;height:4.0pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54572">GSE54572</a><o:p></o:p></span></p>
  </td>
  
 </tr>
 <tr style='mso-yfti-irow:20;height:25.1pt'>
  <td width=51 style='width:50.5pt;border-top:none;border-left:solid black 1.0pt;
  border-bottom:solid windowtext 1.0pt;border-right:solid black 1.0pt;
  mso-border-left-alt:solid black 1.0pt;mso-border-bottom-alt:solid windowtext .5pt;
  mso-border-right-alt:solid black 1.0pt;padding:0in 1.4pt 0in 1.4pt;
  height:25.1pt'>
  <p class=MsoNormal align=center style='text-align:center'><b><span
  style='font-size:8.0pt;color:black'>AAD<o:p></o:p></span></b></p>
  </td>
  <td width=55 style='width:55.25pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid black 1.0pt;
  mso-border-bottom-alt:solid windowtext .5pt;mso-border-right-alt:solid black 1.0pt;
  background:white;padding:0in 1.4pt 0in 1.4pt;height:25.1pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>17<o:p></o:p></span></p>
  </td>
  <td width=54 style='width:.75in;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid black 1.0pt;mso-border-bottom-alt:
  solid windowtext .5pt;mso-border-right-alt:solid black 1.0pt;background:white;
  padding:0in 1.4pt 0in 1.4pt;height:25.1pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>15<o:p></o:p></span></p>
  </td>
  <td width=72 style='width:72.3pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid black 1.0pt;mso-border-bottom-alt:
  solid windowtext .5pt;mso-border-right-alt:solid black 1.0pt;background:white;
  padding:0in 1.4pt 0in 1.4pt;height:25.1pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Superior frontal cortex<o:p></o:p></span></p>
  </td>
  <td width=103 style='width:102.75pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid black 1.0pt;
  mso-border-bottom-alt:solid windowtext .5pt;mso-border-right-alt:solid black 1.0pt;
  background:white;padding:0in 1.4pt 0in 1.4pt;height:25.1pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Illumina HumanHT-12 V3<o:p></o:p></span></p>
  </td>
  <td width=55 style='width:55.2pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid black 1.0pt;mso-border-bottom-alt:
  solid windowtext .5pt;mso-border-right-alt:solid black 1.0pt;background:white;
  padding:0in 1.4pt 0in 1.4pt;height:25.1pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>Mayfield<o:p></o:p></span></p>
  </td>
  <td width=74 style='width:74.0pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid black 1.0pt;mso-border-bottom-alt:
  solid windowtext .5pt;mso-border-right-alt:solid black 1.0pt;background:white;
  padding:0in 1.4pt 0in 1.4pt;height:25.1pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29555">GSE29555</a><o:p></o:p></span></p>
  </td>
  
 </tr>
 <tr style='mso-yfti-irow:21;height:18.85pt'>
  <td width=51 style='width:50.5pt;border-top:none;border-left:solid black 1.0pt;
  border-bottom:solid windowtext 1.0pt;border-right:solid black 1.0pt;
  mso-border-left-alt:solid black 1.0pt;mso-border-bottom-alt:solid windowtext .5pt;
  mso-border-right-alt:solid black 1.0pt;padding:0in 1.4pt 0in 1.4pt;
  height:18.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><b><span
  style='font-size:8.0pt;color:black'>TOTAL<o:p></o:p></span></b></p>
  </td>
  <td width=55 style='width:55.25pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid black 1.0pt;
  mso-border-bottom-alt:solid windowtext .5pt;mso-border-right-alt:solid black 1.0pt;
  background:white;padding:0in 1.4pt 0in 1.4pt;height:18.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>407<o:p></o:p></span></p>
  </td>
  <td width=54 style='width:.75in;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid black 1.0pt;mso-border-bottom-alt:
  solid windowtext .5pt;mso-border-right-alt:solid black 1.0pt;background:white;
  padding:0in 1.4pt 0in 1.4pt;height:18.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'>293<o:p></o:p></span></p>
  </td>
  <td width=72 style='width:72.3pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid black 1.0pt;mso-border-bottom-alt:
  solid windowtext .5pt;mso-border-right-alt:solid black 1.0pt;background:white;
  padding:0in 1.4pt 0in 1.4pt;height:18.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'><o:p>&nbsp;</o:p></span></p>
  </td>
  <td width=103 style='width:102.75pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid black 1.0pt;
  mso-border-bottom-alt:solid windowtext .5pt;mso-border-right-alt:solid black 1.0pt;
  background:white;padding:0in 1.4pt 0in 1.4pt;height:18.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'><o:p>&nbsp;</o:p></span></p>
  </td>
  <td width=55 style='width:55.2pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid black 1.0pt;mso-border-bottom-alt:
  solid windowtext .5pt;mso-border-right-alt:solid black 1.0pt;background:white;
  padding:0in 1.4pt 0in 1.4pt;height:18.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'><o:p>&nbsp;</o:p></span></p>
  </td>
  <td width=74 style='width:74.0pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid black 1.0pt;mso-border-bottom-alt:
  solid windowtext .5pt;mso-border-right-alt:solid black 1.0pt;background:white;
  padding:0in 1.4pt 0in 1.4pt;height:18.85pt'>
  <p class=MsoNormal align=center style='text-align:center'><span
  style='font-size:8.0pt;color:black'><o:p>&nbsp;</o:p></span></p>
  </td></table>



### RNAseq datasets used in this study

| Disease | Dataset | # Samples (Cases/Controls) | Brain Region | Library | Sequencer | Data |
| -----   | ------- | ------------------ | -------------|----------|------|----|
| **SCZ** | PsychEncode-GVEX | 53 / 53 | Frontal Cortex | TruSeq / RiboZero | Illumina HiSeq2000 | <a href="https://www.synapse.org/#!Synapse:syn4590909">syn4590909</a>| 
| **BD** | PsychEncode-GVEX | 47 / 53 (same) | Frontal Cortex | TruSeq / RiboZero | Illumina HiSeq2000 | <a href="https://www.synapse.org/#!Synapse:syn4590909">syn4590909</a>| 
| **ASD** |ASD-pancortical | 53 / 35 | BA4/6; BA38; BA7; BA17 | TruSeq / RiboZero | Illumina HiSeq 2500  | <a href="https://www.synapse.org/#!Synapse:syn11242290">syn11242290</a>|
| **SCZ** | CommonMind | 262 / 295 | DLPFC | TruSeq / RiboZero | Illumina HiSeq 2500  | <a href="https://www.synapse.org/#!Synapse:syn2759792">syn2759792</a>|
| **BD** |CommonMind | 47 / 295 (same) | DLPFC | TruSeq / RiboZero | Illumina HiSeq 2500  | <a href="https://www.synapse.org/#!Synapse:syn2759792">syn2759792</a>|


### GWAS datasets used in this study

<table class=MsoNormalTable border=0 cellspacing=0 cellpadding=0 width=418
 style='border-collapse:collapse;mso-table-layout-alt:fixed;border:none;
 mso-border-alt:solid windowtext .5pt;mso-yfti-tbllook:1184;mso-padding-alt:
 0in 5.4pt 0in 5.4pt;mso-border-insideh:.5pt solid windowtext;mso-border-insidev:
 .5pt solid windowtext'>
 <tr style='mso-yfti-irow:0;mso-yfti-firstrow:yes;height:20.65pt'>
  <td width=63 style='width:62.75pt;border:solid windowtext 1.0pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt;height:20.65pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><b style='mso-bidi-font-weight:normal'><span lang=EN-GB
  style='font-size:10.0pt;mso-fareast-font-family:Calibri;mso-ansi-language:
  EN-GB'>Disorder / Trait<o:p></o:p></span></b></p>
  </td>
  <td width=79 style='width:78.9pt;border:solid windowtext 1.0pt;border-left:
  none;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0in 5.4pt 0in 5.4pt;height:20.65pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><b style='mso-bidi-font-weight:normal'><span lang=EN-GB
  style='font-size:10.0pt;mso-fareast-font-family:Calibri;mso-ansi-language:
  EN-GB'>Consortium<o:p></o:p></span></b></p>
  </td>
  <td width=164 style='width:164.0pt;border:solid windowtext 1.0pt;border-left:
  none;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0in 5.4pt 0in 5.4pt;height:20.65pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><b style='mso-bidi-font-weight:normal'><span lang=EN-GB
  style='font-size:10.0pt;mso-fareast-font-family:Calibri;mso-ansi-language:
  EN-GB'>Dataset<o:p></o:p></span></b></p>
  </td>
  <td width=36 style='width:36.1pt;border:solid windowtext 1.0pt;border-left:
  none;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0in 5.4pt 0in 5.4pt;height:20.65pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><b style='mso-bidi-font-weight:normal'><span lang=EN-GB
  style='font-size:10.0pt;mso-fareast-font-family:Calibri;mso-ansi-language:
  EN-GB'>Date<o:p></o:p></span></b></p>
  </td>
  <td width=77 style='width:76.5pt;border:solid windowtext 1.0pt;border-left:
  none;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0in 5.4pt 0in 5.4pt;height:20.65pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><b style='mso-bidi-font-weight:normal'><span lang=EN-GB
  style='font-size:10.0pt;mso-fareast-font-family:Calibri;mso-ansi-language:
  EN-GB'>Total Sample <br>
  Size (Cases)<o:p></o:p></span></b></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:1;height:20.2pt'>
  <td width=63 style='width:62.75pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0in 5.4pt 0in 5.4pt;height:20.2pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>SCZ<o:p></o:p></span></p>
  </td>
  <td width=79 style='width:78.9pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt;height:20.2pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>PGC<o:p></o:p></span></p>
  </td>
  <td width=164 style='width:164.0pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt;height:20.2pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><a
  href="https://www.med.unc.edu/pgc/files/resultfiles/scz2.snp.results.txt.gz"><span
  lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:Calibri;
  mso-ansi-language:EN-GB'>SCZ2.snp.results.txt.gz</span></a><span lang=EN-GB
  style='font-size:10.0pt;mso-fareast-font-family:Calibri;mso-ansi-language:
  EN-GB'><o:p></o:p></span></p>
  </td>
  <td width=36 style='width:36.1pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt;height:20.2pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>2014<o:p></o:p></span></p>
  </td>
  <td width=77 style='width:76.5pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt;height:20.2pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>82,315<br>
  (35,476)<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:2'>
  <td width=63 style='width:62.75pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>BD<o:p></o:p></span></p>
  </td>
  <td width=79 style='width:78.9pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>PGC<o:p></o:p></span></p>
  </td>
  <td width=164 style='width:164.0pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt'>
  <p class=MsoNormal align=center style='text-align:center'><a
  href="https://www.med.unc.edu/pgc/files/resultfiles/pgc.bip.2012-04.zip"><span
  style='font-size:10.0pt;color:#994100;background:white'>pgc.bip.2012-04.zip</span></a><span
  style='font-size:10.0pt'><o:p></o:p></span></p>
  </td>
  <td width=36 style='width:36.1pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>2012<o:p></o:p></span></p>
  </td>
  <td width=77 style='width:76.5pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>16,731<br>
  (7,481)<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:3'>
  <td width=63 style='width:62.75pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>MDD<o:p></o:p></span></p>
  </td>
  <td width=79 style='width:78.9pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>PGC<o:p></o:p></span></p>
  </td>
  <td width=164 style='width:164.0pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><a
  href="https://www.med.unc.edu/pgc/files/resultfiles/pgc.mdd.2012-04.zip"><span
  lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:Calibri;
  mso-ansi-language:EN-GB'>pgc.mdd.2012-04.zip</span></a><span lang=EN-GB
  style='font-size:10.0pt;mso-fareast-font-family:Calibri;mso-ansi-language:
  EN-GB'><o:p></o:p></span></p>
  </td>
  <td width=36 style='width:36.1pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>2012<o:p></o:p></span></p>
  </td>
  <td width=77 style='width:76.5pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>18,759<br>
  (9,240)<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:4'>
  <td width=63 style='width:62.75pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>ASD<o:p></o:p></span></p>
  </td>
  <td width=79 style='width:78.9pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span class=SpellE><span lang=EN-GB style='font-size:10.0pt;
  mso-fareast-font-family:Calibri;mso-ansi-language:EN-GB'>iPSYCH</span></span><span
  lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:Calibri;
  mso-ansi-language:EN-GB'><o:p></o:p></span></p>
  </td>
  <td width=164 style='width:164.0pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span style='font-size:10.0pt;mso-fareast-font-family:Calibri'><o:p><a href="./raw_data/GWAS">Summary Statistics</a></o:p></span></p>
  </td>
  <td width=36 style='width:36.1pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>2017<o:p></o:p></span></p>
  </td>
  <td width=77 style='width:76.5pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>28,131<br>
  (8,605)<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:5'>
  <td width=63 style='width:62.75pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>ASD<o:p></o:p></span></p>
  </td>
  <td width=79 style='width:78.9pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>PGC<o:p></o:p></span></p>
  </td>
  <td width=164 style='width:164.0pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><a
  href="https://www.med.unc.edu/pgc/files/resultfiles/pgcasdeuro.gz"><span
  style='font-size:10.0pt;mso-fareast-font-family:Calibri;background:white'>PGC.ASD.euro.all.25Mar2015.txt.gz</span></a><span
  lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:Calibri;
  mso-ansi-language:EN-GB'><o:p></o:p></span></p>
  </td>
  <td width=36 style='width:36.1pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>2015<o:p></o:p></span></p>
  </td>
  <td width=77 style='width:76.5pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>10,610<br>
  (5305)<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:6'>
  <td width=63 style='width:62.75pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>AAD<o:p></o:p></span></p>
  </td>
  <td width=79 style='width:78.9pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span class=SpellE><span lang=EN-GB style='font-size:10.0pt;
  mso-fareast-font-family:Calibri;mso-ansi-language:EN-GB'>AlcGen</span></span><span
  lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:Calibri;
  mso-ansi-language:EN-GB'><o:p></o:p></span></p>
  </td>
  <td width=164 style='width:164.0pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'><o:p><a href="https://www.ncbi.nlm.nih.gov/pubmed/21471458">Schumann et al., 2011</a></o:p></span></p>
  </td>
  <td width=36 style='width:36.1pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>2011<o:p></o:p></span></p>
  </td>
  <td width=77 style='width:76.5pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;color:black;mso-ansi-language:EN-GB'>23,347 <br>
  (NA)<span style='background:yellow;mso-highlight:yellow'><o:p></o:p></span></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:7;height:17.5pt'>
  <td width=63 style='width:62.75pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0in 5.4pt 0in 5.4pt;height:17.5pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>IBD<o:p></o:p></span></p>
  </td>
  <td width=79 style='width:78.9pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt;height:17.5pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>IIBDGC<o:p></o:p></span></p>
  </td>
  <td width=164 style='width:164.0pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt;height:17.5pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span style='font-size:10.0pt;mso-fareast-font-family:Calibri'>EUR.IBD.gwas.assoc.txt</span><span
  lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:Calibri;
  mso-ansi-language:EN-GB'><o:p></o:p></span></p>
  </td>
  <td width=36 style='width:36.1pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt;height:17.5pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>2015<o:p></o:p></span></p>
  </td>
  <td width=77 style='width:76.5pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt;height:17.5pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>34,652<br>
  (12,882)<span style='background:yellow;mso-highlight:yellow'><o:p></o:p></span></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:8;height:17.5pt'>
  <td width=63 style='width:62.75pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0in 5.4pt 0in 5.4pt;height:17.5pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>Educational Attainment<o:p></o:p></span></p>
  </td>
  <td width=79 style='width:78.9pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt;height:17.5pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>SSGAC<o:p></o:p></span></p>
  </td>
  <td width=164 style='width:164.0pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt;height:17.5pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><a href="http://ssgac.org/documents/EduYears_Main.txt.gz"><span
  style='font-size:10.0pt;mso-fareast-font-family:Calibri'>EduYears_Main.txt.gz</span></a><span
  style='font-size:10.0pt;mso-fareast-font-family:Calibri'><o:p></o:p></span></p>
  </td>
  <td width=36 style='width:36.1pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt;height:17.5pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>2016<o:p></o:p></span></p>
  </td>
  <td width=77 style='width:76.5pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt;height:17.5pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>328,917<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:9;height:16.6pt'>
  <td width=63 style='width:62.75pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0in 5.4pt 0in 5.4pt;height:16.6pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>Subjective Well Being<o:p></o:p></span></p>
  </td>
  <td width=79 rowspan=3 style='width:78.9pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt;height:16.6pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>SSGAC<o:p></o:p></span></p>
  </td>
  <td width=164 style='width:164.0pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt;height:16.6pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><a href="http://ssgac.org/documents/SWB_Full.txt.gz"><span
  style='font-size:10.0pt;mso-fareast-font-family:Calibri'>SWB_Full.txt.gz</span></a><span
  style='font-size:10.0pt;mso-fareast-font-family:Calibri'><o:p></o:p></span></p>
  </td>
  <td width=36 rowspan=3 style='width:36.1pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt;height:16.6pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>2016<o:p></o:p></span></p>
  </td>
  <td width=77 style='width:76.5pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt;height:16.6pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>298,420<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:10;height:16.6pt'>
  <td width=63 style='width:62.75pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0in 5.4pt 0in 5.4pt;height:16.6pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>Neuroticism<o:p></o:p></span></p>
  </td>
  <td width=164 style='width:164.0pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt;height:16.6pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><a href="http://ssgac.org/documents/Neuroticism_Full.txt.gz"><span
  style='font-size:10.0pt;mso-fareast-font-family:Calibri'>Neuroticism_Full.txt.gz</span></a><span
  style='font-size:10.0pt;mso-fareast-font-family:Calibri'><o:p></o:p></span></p>
  </td>
  <td width=77 style='width:76.5pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt;height:16.6pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>170,911<o:p></o:p></span></p>
  </td>
 </tr>
 <tr style='mso-yfti-irow:11;mso-yfti-lastrow:yes;height:16.6pt'>
  <td width=63 style='width:62.75pt;border:solid windowtext 1.0pt;border-top:
  none;mso-border-top-alt:solid windowtext .5pt;mso-border-alt:solid windowtext .5pt;
  padding:0in 5.4pt 0in 5.4pt;height:16.6pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>Depressive Symptoms<o:p></o:p></span></p>
  </td>
  <td width=164 style='width:164.0pt;border-top:none;border-left:none;
  border-bottom:solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;
  mso-border-top-alt:solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;
  mso-border-alt:solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt;height:16.6pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><a href="http://ssgac.org/documents/DS_Full.txt.gz"><span
  style='font-size:10.0pt;mso-fareast-font-family:Calibri'>DS_Full.txt.gz</span></a><span
  style='font-size:10.0pt;mso-fareast-font-family:Calibri'><o:p></o:p></span></p>
  </td>
  <td width=77 style='width:76.5pt;border-top:none;border-left:none;border-bottom:
  solid windowtext 1.0pt;border-right:solid windowtext 1.0pt;mso-border-top-alt:
  solid windowtext .5pt;mso-border-left-alt:solid windowtext .5pt;mso-border-alt:
  solid windowtext .5pt;padding:0in 5.4pt 0in 5.4pt;height:16.6pt'>
  <p class=Paragraph align=center style='margin-top:0in;text-align:center;
  text-indent:0in'><span lang=EN-GB style='font-size:10.0pt;mso-fareast-font-family:
  Calibri;mso-ansi-language:EN-GB'>161,460<o:p></o:p></span></p>
  </td>
 </tr>
</table>
