void tmp_program()
{
TRandom3 *random = new TRandom3();

std::cout<<"seed of random number before "<<random->GetSeed()<<std::endl;
random->SetSeed(5);

std::cout<<"seed of random number"<<random->GetSeed()<<std::endl;

float random_number;

for(int i=0; i<100; i++){
//random_number = random->Gaus(0,1);
random_number = random->Rndm();
std::cout<<random_number<<std::endl;
}
}

