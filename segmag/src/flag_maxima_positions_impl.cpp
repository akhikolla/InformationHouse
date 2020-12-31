#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector flag_maxima_positions_impl(NumericVector values)
{  
  int n_values = values.size();
  NumericVector maxima(n_values);

  /** Walks through the array "values" and flags local maxima in values with a
   *  1 in the array "maxima"
   * 
   *  If multiple values have the same value at the maximum, then the center of
   *  the maximum is defined as the center of these values. If the point of the
   *  maximum is a .5 value than it is rounded up (e.g. 6,7,7,4 --> 0,0,1,0)
   */

  for (int i = 0; i < n_values; i++)
  {
    // Sonderbehandlung linke Seite
    if (i == 0)
    {
      // vorwärts suchen, wann erster Wert != values[0] findet
      for (int j = 1; j < n_values; j++)
      {
        if (values[0] > values[j])
        {
          // Ja, war Maximum da nach rechts zu fallen beginnt
          // Die Mitte der Spitze als Maximum markieren
          maxima[std::ceil((j-1)/2.0)] = 1;
          break;
        }
        else if (values[0] < values[j])
        {
          // Nein, Funktion steigt nach rechts weiter, war kein Maximum an i==0 --> Suche abbrechen
          break;
        }
      }
    }
    else
    {
      // Sonderbehandlung rechte Seite
      if (i == n_values - 1)
      {
        // Zurück suchen, wann erster Wert != values[i] findet
        for (int j = i-1; j >= 0; j--)
        {
          if (values[i] > values[j])
          {
            // Ja, war Maximum da wieder zu fallen beginnt
            // Die Mitte der Spitze als Maximum markieren
            maxima[std::ceil((i+j+1)/2.0)] = 1;
            break;
          }
          else if (values[i] < values[j])
          {
            // Nein, Funktion steigt nach links weiter, war kein Maximum an i --> Suche abbrechen
            break;
          }
        }
      }
      
      // Normal: Potentielles Maximum bei i-1: Werte beginnen zu fallen;
      // ACHTUNG1: kein else if, da sonst folgender Fall nicht funktioniert c(1,2,1), Grund: hier wird i-1 geprüft und oben wird i (n_values-1) geprüft
      // ACHTUNG2: hier in else branch drin, da i > 0 Voraussetzung für Prüfung values[i-1] ist      
      if (values[i] < values[i-1])
      {      
        // Zurück suchen, wann erster Wert != values[i-1] findet
        for (int j = i-2; j >= 0; j--)
        {
          if (values[i-1] > values[j])
          {
            // Ja, war Maximum da wieder zu fallen beginnt
            // Die Mitte der Spitze als Maximum markieren
            // ((i-1) + (j+1)) / 2
            maxima[std::ceil((i+j)/2.0)] = 1;
            break;
          }
          else if (values[i-1] < values[j])
          {
            // Nein, Funktion steigt nach links weiter, war kein Maximum an i-1 --> Suche abbrechen
            break;
          }
        }
      }
    }
  }
  
  return( maxima );
}
