{
  description = "Robust Implementation of the LaBRADOR Proof System";

  inputs = {
    nixpkgs.url = "https://flakehub.com/f/NixOS/nixpkgs/*.tar.gz";
    cargo2nix.url = "github:cargo2nix/cargo2nix/release-0.11.0";
    flake-utils.url = "github:numtide/flake-utils/v1.0.0";
  };

  outputs = inputs: with inputs;
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          overlays = [ cargo2nix.overlays.default ];
        };
        inherit (pkgs) lib;

        rustPackageSet = pkgs.rustBuilder.makePackageSet {
          rustVersion = "1.71.1";
          extraRustComponents = [ "rustfmt" "clippy" ];
        };

        buildInputs = [
          pkgs.cargo-all-features
          pkgs.cargo-deny
          pkgs.cargo-nextest
          pkgs.rustup

        workspaceShell = rustPackageSet.workspaceShell {
          packages = buildInputs;
        };
      in
      rec
      {
        packages = {
          default = mps { };
          tests = mps { compileMode = "test"; };
          ci = pkgs.rustBuilder.runTests mps {
            RUST_BACKTRACE = "full";
          };
        };

        devShell = workspaceShell;
      }
    );
}
