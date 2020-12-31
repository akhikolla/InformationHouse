#include <Rcpp.h>
#include "zip.h"
#include "rapidxml.h"
#include "xlsxnames.h"

using namespace Rcpp;

xlsxnames::xlsxnames(const std::string& path) {
  // Names are stored at the workbook level, even if scoped to sheets
  std::string book = zip_buffer(path, "xl/workbook.xml");

  rapidxml::xml_document<> xml;
  xml.parse<rapidxml::parse_strip_xml_namespaces>(&book[0]);

  rapidxml::xml_node<>* workbook = xml.first_node("workbook");
  rapidxml::xml_node<>* definedNames = workbook->first_node("definedNames");

  int n(0);
  if (definedNames != NULL) {
    // Count the rules.
    for (rapidxml::xml_node<>* definedName = definedNames->first_node("definedName");
        definedName; definedName = definedName->next_sibling("definedName")) {
        ++n;
    }
  }

  name_     = CharacterVector(n, NA_STRING);
  sheet_id_ = IntegerVector(n,   NA_INTEGER);
  formula_  = CharacterVector(n, NA_STRING);
  comment_  = CharacterVector(n, NA_STRING);
  hidden_   = LogicalVector(n,   NA_LOGICAL);

  if (definedNames != NULL) {
    int i(0);
    for (rapidxml::xml_node<>* definedName = definedNames->first_node("definedName");
        definedName; definedName = definedName->next_sibling("definedName")) {

      rapidxml::xml_attribute<>* name = definedName->first_attribute("name");
      if (name != NULL) {
        name_[i] = name->value();
      }

      rapidxml::xml_attribute<>* localSheetId = definedName->first_attribute("localSheetId");
      if (localSheetId != NULL) {
        sheet_id_[i] = strtol(localSheetId->value(), NULL, 10);
      }

      formula_[i] = definedName->value();

      rapidxml::xml_attribute<>* comment = definedName->first_attribute("comment");
      if (comment != NULL) {
        comment_[i] = comment->value();
      }

      rapidxml::xml_attribute<>* hidden = definedName->first_attribute("hidden");
      if (hidden == NULL) {
        hidden_[i] = false;
      } else { // no known case
        std::string hidden_value(hidden->value()); // # nocov start
        if (hidden_value == "1" || hidden_value == "true") {
          hidden_[i] = true;
        } else {
          hidden_[i] = false;
        } // # nocov end
      }

      ++i;
    }
  }
}

List& xlsxnames::information() {
  // Returns a nested data frame of everything, the data frame itself wrapped in
  // a list.
  information_ = List::create(
      _["sheet_id"] = sheet_id_,
      _["name"] =     name_,
      _["formula"] =  formula_,
      _["comment"] =  comment_,
      _["hidden"] =   hidden_);

  // Turn list of vectors into a data frame without checking anything
  int n = Rf_length(information_[0]);
  information_.attr("class") = CharacterVector::create("tbl_df", "tbl", "data.frame");
  information_.attr("row.names") = IntegerVector::create(NA_INTEGER, -n); // Dunno how this works (the -n part)

  return information_;
}
