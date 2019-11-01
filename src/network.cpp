#include "network.h"



//Calculates the number and total intensity of connections to neuron \p n.
  //\param n : the index of the receiving neuron. n==link.first.second
  //\return a pair {number of connections, sum of link intensities}.
  // appel: link.first , link.second, link.first.first ect...
 //link = std::pair<const std::pair<long unsigned int, long unsigned int> double>; // pair2=(premier int # neuron1, 2eme int # neuron2) , pair1=(pair2, intensité)


std::pair<size_t, double> Network::degree(const size_t& n) const
{ 
	 	size_t nbr_connections(neighbors(n).size());
		double intensity_sum(0.0);
	 
	 for (auto& links : neighbors(n)){
		
			 intensity_sum += links.second;
		 
	 } 
 std::pair<size_t, double> paire (nbr_connections, intensity_sum);
 return paire;
 }

//firing je le reset mais stock les valeur-etat
// distinction entre inhib et exitateur et voir si ils sont firing avant de sommer leur intensité. Pour les neurones voisins.



 //Performs one time-step of the simulation.
  //\param input : a vector of random values as thalamic input, one value for each neuron. The variance of these values corresponds to excitatory neurons.
  //\return the indices of firing neurons
  //variance de vect[i] = w du pdf
  
std::set<size_t> Network::step(const std::vector<double>& thalam)
{  
  std::set<size_t> firing;
  std::vector<bool> firing_stock(neurons.size());
  
  double I, sum_inhib,sum_ex;
  
  for(size_t i(0); i < neurons.size(); i++){ //itere sur le vect de neurones
	  
	  if(neurons[i].firing()){
		  firing.insert(i);
		  firing_stock.push_back(true); //si il est firing je le stock et le reset
		  neurons[i].reset(); 
			
		}
		else{
			firing_stock.push_back(false);
		}
		
		
		for (auto& links : neighbors(i)){
			
				if(firing_stock[links.first]){
					
					if (neurons[links.first].is_inhibitory ())	// on ajoute soutrait de intensity_sum 
					{
						sum_inhib += links.second;
					}
					else { 													// on ajoute 0.5 intensity_sum de I
						sum_ex += links.second;
					} 
				}
		} 
		
		
		  if (neurons[i].is_inhibitory()){
			  
			I = 0.4*thalam[i] - sum_inhib + 0.5*sum_ex; //on ajoute parce que askip c'est deja negatif
			
		}else {
			I = thalam[i] - sum_inhib + 0.5*sum_ex;
			}
		
		neurons[i].input(I);
		neurons[i].step();
		
	}
 
  return firing;
}


 //Finds the list of neurons with incoming connections to \p n.
  //\param n : the index of the receiving neuron. n==link.first.second
  //\return a vector of pairs {neuron index, link intensity}

std::vector<std::pair<size_t, double> > Network::neighbors(const size_t& n) const
{
  std::pair<size_t, double> index_intens;
  std::vector<std::pair<size_t, double> > neighb;
  
  for (auto i: links){
	  if (n == i.first.first){
		 
		  index_intens.first = i.first.second;
		  index_intens.second = i.second;
		  neighb.push_back(index_intens);
	  }  
  }
  
  return neighb;
}

void Network::resize(const size_t &n, double inhib) {
    size_t old = size();
    neurons.resize(n);
    if (n <= old) return;
    size_t nfs(inhib*(n-old)+.5);
    set_default_params({{"FS", nfs}}, old);
}

void Network::set_default_params(const std::map<std::string, size_t> &types,
                                 const size_t start) {
    size_t k(0), ssize(size()-start), kmax(0);
    std::vector<double> noise(ssize);
    _RNG->uniform_double(noise);
    for (auto I : types) 
        if (Neuron::type_exists(I.first)) 
            for (kmax+=I.second; k<kmax && k<ssize; k++) 
                neurons[start+k].set_default_params(I.first, noise[k]);
    for (; k<ssize; k++) neurons[start+k].set_default_params("RS", noise[k]);
}

void Network::set_types_params(const std::vector<std::string> &_types,
                               const std::vector<NeuronParams> &_par,
                               const size_t start) {
    for (size_t k=0; k<_par.size(); k++) {
        neurons[start+k].set_type(_types[k]);
        neurons[start+k].set_params(_par[k]);
    }
}

void Network::set_values(const std::vector<double> &_poten, const size_t start) {
    for (size_t k=0; k<_poten.size(); k++) 
        neurons[start+k].potential(_poten[k]);
}

bool Network::add_link(const size_t &a, const size_t &b, double str) {
    if (a==b || a>=size() || b>=size() || str<1e-6) return false;
    if (links.count({a,b})) return false;
    if (neurons[b].is_inhibitory()) str *= -2.0;
    links.insert({{a,b}, str});
    return true;
}

size_t Network::random_connect(const double &mean_deg, const double &mean_streng) {
    links.clear();
    std::vector<int> degrees(size());
    _RNG->poisson(degrees, mean_deg);
    size_t num_links = 0;
    std::vector<size_t> nodeidx(size());
    std::iota(nodeidx.begin(), nodeidx.end(), 0);
    for (size_t node=0; node<size(); node++) {
        _RNG->shuffle(nodeidx);
        std::vector<double> strength(degrees[node]);
        _RNG->uniform_double(strength, 1e-6, 2*mean_streng);
        int nl = 0;
        for (size_t nn=0; nn<size() && nl<degrees[node]; nn++)
            if (add_link(node, nodeidx[nn], strength[nl])) nl++;
        num_links += nl;
    }
    return num_links;
}

std::vector<double> Network::potentials() const {
    std::vector<double> vals;
    for (size_t nn=0; nn<size(); nn++)
        vals.push_back(neurons[nn].potential());
    return vals;
}

std::vector<double> Network::recoveries() const {
    std::vector<double> vals;
    for (size_t nn=0; nn<size(); nn++)
        vals.push_back(neurons[nn].recovery());
    return vals;
}

void Network::print_params(std::ostream *_out) {
    (*_out) << "Type\ta\tb\tc\td\tInhibitory\tdegree\tvalence" << std::endl;
    for (size_t nn=0; nn<size(); nn++) {
        std::pair<size_t, double> dI = degree(nn);
        (*_out) << neurons[nn].formatted_params() 
                << '\t' << dI.first << '\t' << dI.second
                << std::endl;
    }
}

void Network::print_head(const std::map<std::string, size_t> &_nt, 
                         std::ostream *_out) {
    size_t total = 0;
    for (auto It : _nt) {
        total += It.second;
        for (auto In : neurons)
            if (In.is_type(It.first)) {
                (*_out) << '\t' << It.first << ".v"
                        << '\t' << It.first << ".u"
                        << '\t' << It.first << ".I";
                break;
            }
    }
    if (total<size())
        for (auto In : neurons) 
            if (In.is_type("RS")) {
                (*_out) << '\t' << "RS.v" << '\t' << "RS.u" << '\t' << "RS.I";
                break;
            }
    (*_out) << std::endl;
}

void Network::print_traj(const int time, const std::map<std::string, size_t> &_nt, 
                         std::ostream *_out) {
    (*_out)  << time;
    size_t total = 0;
    for (auto It : _nt) {
        total += It.second;
        for (auto In : neurons) 
            if (In.is_type(It.first)) {
                (*_out) << '\t' << In.formatted_values();
                break;
            }
    }
    if (total<size())
        for (auto In : neurons) 
            if (In.is_type("RS")) {
                (*_out) << '\t' << In.formatted_values();
                break;
            }
    (*_out) << std::endl;
}
