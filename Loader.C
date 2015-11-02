{
  gROOT->ProcessLine(".L HH4bAna_v2.C+");
  HH4bAna_v2 t;
  t.GetEntry(12); // Fill t data members with entry number 12
  t.Show();       // Show values of entry 12
  t.Show(16);     // Read and show values of entry 16
  t.Loop();

}
