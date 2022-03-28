<script type="text/javascript">
  // When the document is fully rendered...
  $(document).ready(function() {
    // ...select all header elements...
    $('h1, h2, h3, h4, h5').each(function() {
      // ...and add an id to them corresponding to their 'titles'
      $(this).attr('id', $(this).html());
    });
  });
</script>
# TCR repertoire analysis

<a href="#TCR analysis section">TCR analysis section</a><br>
<a href="#Motif analysis section">Motif analysis section</a><br>
<a href="#Diversity and chain interrogation">Diversity and chain interrogation</a><br>
<a href="#Group overlap analysis">Group overlap analysis</a>

## Side panel.
Upload the file. This can be from our QC section or alternative sources. 
The other features in the side panel are 
- 'Type of group'
  + This is used to change the comparison. We recommend either using "group","indiv" or "group.indiv". 
- 'Type of data'
  + This segregates out if the original file was 'raw' or 'summarized'
- 'Type of font' 
  + Specify the font for the figures. the R default fonts are serif, sans and mono. Additional fonts were found on https://fonts.google.com (email Kerry if there is a specific font you would like to use.)

<img src="inst/extdata/Images/upload.TCR.analysis.png" width="300">
<img src="inst/extdata/Images/side.panel.png" width="300">

## TCR analysis section
<a href="#summary table">summary table</a><br>
<a href="#Treemap">Treemap</a><br>
<a href="#Chord diargram">Chord diargram</a><br>
<a href="#Pie chart">Pie chart</a><br>

### Overview of TCR pairing

#### summary table

The user can specify the type of summary table to download. 

They can either select their own columns (general summary) or downlaod as TCRdist3 .csv output.

<img src="inst/extdata/Images/general.summary.png" width="600">

For the TCRdist3, there is a need to use our QC process as it matches the IMGT column names. 

There is also a need to select if the input data is either alpha-beta (ab) or gamma-delta (gd) for the TCRdist3 column selection.

<img src="inst/extdata/Images/TCRdist3.summary.png" width="600">

<a href="#TCR analysis section">TCR analysis section</a><br>
<a href="#TCR repertoire analysis">Go to top</a><br>

#### Treemap
The user can specify: 
- The order of the group (i.e. CD8 and IFNg)
- colour choices include: default, rainbow, random or one colour (specified in side panel) 
    + The colour can be altered afterwards
- If they want the labels to appear on the graph
- Column to colour as well as column to separate the panel
- This plot can be downloaded as a PNG or PDF

<img src="inst/extdata/Images/treemap.png" width="800">

<a href="#TCR analysis section">TCR analysis section</a><br>
<a href="#TCR repertoire analysis">Go to top</a><br>

#### Chord diargram

There are several features the user can specify:
- Sub-group to display
- Which 
- Colour choices: default, rainbow, random or one colour (specified in side panel) 
- 


<img src="inst/extdata/Images/chord.png" width="800">

<a href="#TCR analysis section">TCR analysis section</a><br>
<a href="#TCR repertoire analysis">Go to top</a><br>

#### Pie chart




<img src="inst/extdata/Images/pie.png" width="800">

<a href="#TCR analysis section">TCR analysis section</a><br>
<a href="#TCR repertoire analysis">Go to top</a><br>

### Motif analysis section

<a href="#TCR analysis section">TCR analysis section</a><br>
<a href="#TCR repertoire analysis">Go to top</a><br>

#### 

the nucleotide and amino acid plots show the unique sequences of a certain length

<img src="IMAGES/aaMotif.png" width="1000">

### Diversity and chain interrogation

<a href="#TCR analysis section">TCR analysis section</a><br>
<a href="#TCR repertiore analysis">Go to top</a><br>

### Group overlap analysis

<a href="#TCR analysis section">TCR analysis section</a><br>
<a href="#TCR repertiore analysis">Go to top</a><br>
