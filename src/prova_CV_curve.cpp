/// CV curve (std::vector<double>& phi0, double Vgstart, double Vgend, double dVg)
// Vgate da -30 a 30 [V] con dVgate = 0.1

double	dVg = 0.1,
		Vgstart = -30,
		Vgend = 30,
		Vg = 0;
std::vector<double> Vguess = phi0;
		
for(auto i=Vgstart; i<=Vgend; i+=dVg){
	_VG = i;
	P.NLPoisson(Vguess,true);	/// in realtà ha condizioni al contorno differenti!!!
	
	Vguess = P._Vin;
}

risolvo NLPoisson con Vguess phi0 e Vg fissata.
vario Vgate e uso ogni volta la Vguess calcolata precedentemente