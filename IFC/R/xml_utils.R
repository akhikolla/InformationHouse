################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2020 Yohann Demont                                             #
#                                                                              #
# It is part of IFC package, please cite:                                      #
# -IFC: An R Package for Imaging Flow Cytometry                                #
# -YEAR: 2020                                                                  #
# -COPYRIGHT HOLDERS: Yohann Demont, Gautier Stoll, Guido Kroemer,             #
#                     Jean-Pierre Marolleau, Loïc Garçon,                      #
#                     INSERM, UPD, CHU Amiens                                  #
#                                                                              #
# DISCLAIMER:                                                                  #
# -You are using this package on your own risk!                                #
# -We do not guarantee privacy nor confidentiality.                            #
# -This program is distributed in the hope that it will be useful, but WITHOUT #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        #
# FITNESS FOR A PARTICULAR PURPOSE. In no event shall the copyright holders or #
# contributors be liable for any direct, indirect, incidental, special,        #
# exemplary, or consequential damages (including, but not limited to,          #
# procurement of substitute goods or services; loss of use, data, or profits;  #
# or business interruption) however caused and on any theory of liability,     #
# whether in contract, strict liability, or tort (including negligence or      #
# otherwise) arising in any way out of the use of this software, even if       #
# advised of the possibility of such damage.                                   #
#                                                                              #
# You should have received a copy of the GNU General Public License            #
# along with IFC. If not, see <http://www.gnu.org/licenses/>.                  #
################################################################################

#' @title XML Node List Expansion
#' @description
#' Helper to stringify XML node.
#' @param x value return by xml2::as_list. Default is 5.
#' @param max maximum number of recurrence into subnodes.
#' @keywords internal
expand_list = function(x, max = 5) {
  max = max - 1
  N = names(x)
  ans = lapply(x, FUN=function(n) {
    foo = attributes(n)
    bar = (names(foo) == "names")
    if((length(n) == 0) | (max < 0)) return(foo)
    c(foo[!bar], expand_list(n[bar], max = max))
  })
  names(ans) <- N
  return(ans)
}

#' @title XML Entities Protection
#' @description
#' Helper to escape xml entities.
#' @param text value return by xml2::as_list. Default is 5.
#' @details entities will be replaced by:\cr
#' -& to "&amp;"\cr
#' -> to "&gt;"\cr
#' -< to "&lt;"\cr
#' -' to "&apos;"\cr
#' -" to "&quot;"\cr
#' @return a character vector where xml entities have been escaped.
#' @keywords internal
escape_entities = function(text) {
  text = gsub("&", "&amp;", text)
  text = gsub(">", "&gt;", text)
  text = gsub("<", "&lt;", text)
  text = gsub("'", "&apos;", text)
  gsub('"', '&quot;', text)
}

#' @title XML Node to List Conversion
#' @description
#' Helper to convert xml node to R list.
#' @param x A document, node, or node set.
#' @param max maximum number of recurrence into subnodes. Default is 5.
#' @details it acts as_list but value returned is different, with attributes
#' expanded to sublists rather than recovered as attributes
#' @keywords internal
to_list_node = function(x, max = 5) {
  expand_list(xml2::as_list(x), max = max)
}

#' @title List to XML Node Conversion
#' @description
#' Helper to convert R list to xml node (character representation).
#' @param x a list to convert
#' @param name name of the node to create
#' @param kids a list containing children xml nodes elements (each elements should come from to_xml_list)
#' @param indent indent used for kids when provided. Default is "  ".
#' @param escape escape used for kids when provided. Default is "\\n". 
#' @details it acts as_list but value returned is different, with attributes
#' expanded to sublists rather than recovered as attributes
#' @keywords internal
to_xml_list = function(x, name, kids, indent = "  ", escape = "\n") {
  name = as.character(name); assert(name, len = 1, typ = "character")
  N = names(x)
  if(length(N)==0 & missing(kids)) {
    return(sprintf("<%s>%s</%s>", name, x, name))
  }
  node = do.call(paste0, args = c(list(collapse = " ", lapply(1:length(x), FUN=function(i) sprintf('%s="%s"',N[i],x[[i]])))))
  if(missing(kids)) {
    return(sprintf("<%s %s />", name, node))
  } else {
    K = lapply(kids, FUN =function(k) paste0(indent, k, escape))
    return(sprintf("<%s %s>%s%s</%s>", name, node, escape, paste0(K, collapse=""), name))
  }
}

#' @title List to XML Conversion
#' @description
#' Helper to convert R list to xml node (character representation).
#' @param name name of the node to create.
#' @param attrs a named list of name-value pairs to be used as attributes for the XML node.
#' @param .children a list containing XML node elements or content.
#' @param text the text content for the new XML node.
#' @return an R object that points to the C-level structure instance.
#' @keywords internal
xml_new_node <- function(name, attrs, .children, text, ...) {
  tmp <- xml_new_root("foo")
  if(!missing(text) & !missing(attrs)) stop("should be either a 'text' node or a 'attr' node")
  tmp %>% xml_add_child(.value = name)
  node <- xml_find_first(tmp, xpath = name)
  if(!missing(text)) {
    node %>% xml_set_text(value = text) 
  }
  if(!missing(attrs)) {
    node %>% xml_set_attrs(value = attrs)
  }
  if(!missing(.children)) {
    if(any(class(.children) == "xml_node")) {
      node %>% xml_add_child(.value = .children)
    } else {
      L = length(.children)
      first <- node %>% xml_add_child(.value = "")
        for(i in 1:L) {
          first %>% xml_add_sibling(.value = .children[[i]], .where = "before")
        }
      xml_remove(first)
    }
  }
  return(node)
}
