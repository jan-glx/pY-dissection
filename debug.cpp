
#include <string>
#include <iostream>

extern "C" {
#include "anchor.c"
}

// [[Rcpp::export]]
std::string anchor(std::string seq) {
	char* seq_name = "BLA";

	anchor_(seq_name, seq.c_str(), seq.length(), 1, 0);
	return out;
}

int main(){
	std::cout << anchor("MAFKDTGKTPVEPEVAIHRIRITLTSRNVKSLEKVCADLIRGAKEKNLKVKGPVRMPTKTLRITTRKTPCGEGSKTWDRFQMRIHKRLIDLHSPSEIVKQITSISIEPGVEVEVTIADA");
	std::cout << anchor("MAFKDTGKTPVEPEVAIHRIRITLTSRNVKSLEKVCADLIRGAKEKNLKVKGPVRMPTKTLRITTRKTPCGEGSKTWDRFQMRIHKRLIDLHSPSEIVKQITSISIEPGVEVEVTIADA");
	std::cout << "done";
	std::cin.get();
	return 1;
}