#' Update IGV Sessions
#' 
#' @description 
#' This function allows updating the path to the reference genome and annotations 
#' for each IGV session (in XML format) within a user-specified directory. XML 
#' files are searched for in a recursive manner.
#' 
#' @param directory A `String`. Path to the directory containing the xml files. 
#' @param genome A `String`. Path to the genome on which reads have been aligned.
#' @param annotations A `String`. Path to the annotation file (in GFF or GTF
#' format).
#' 
#' @export
setPaths <- function(
  directory,
  genome,
  annotations) {
  
  genome       <- normalizePath(genome)
  annotations  <- normalizePath(annotations)

  all_xml <- list.files(
    directory,
    recursive = TRUE, 
    pattern = ".xml$"
  )
  
  for (i_xml in seq_along(all_xml)){
    
    igv_session <- xml2::read_xml(all_xml[[i_xml]]) 
    xml2::xml_set_attr(igv_session, attr = "genome", value = genome)
    
    if (xml2::xml_attr(igv_session, attr = "stranded") == "TRUE") {
      xml2::xml_set_attr(xml2::xml_children(xml2::xml_children(igv_session)[1])[10], attr = "path", value = annotations)
      xml2::xml_set_attr(xml2::xml_children(xml2::xml_children(igv_session))[15], att = "id", value = genome)
      xml2::xml_set_attr(xml2::xml_children(xml2::xml_children(igv_session))[16], att = "id", value = annotations)
    } else {
      xml2::xml_set_attr(xml2::xml_children(xml2::xml_children(igv_session)[1])[5], attr = "path", value = annotations)
      xml2::xml_set_attr(xml2::xml_children(xml2::xml_children(igv_session))[10], att = "id", value = genome)
      xml2::xml_set_attr(xml2::xml_children(xml2::xml_children(igv_session))[11], att = "id", value = annotations)
    }
    
    xml2::write_xml(igv_session, all_xml[[i_xml]])
  }
}