fvScalarMatrix phiEqn
(
    fvm::ddt(phi)
  + fvm::div(phi, U)
  - fvm::laplacian(D, phi)
);

phiEqn.solve();