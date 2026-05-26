# Memristor SPICE Simulator (C Implementation)

A custom SPICE-like circuit simulation engine written in C, developed as part of a Diploma Thesis in Electrical and Computer Engineering.

The simulator supports standard circuit analyses and introduces **native memristor modeling**, enabling efficient simulation of **Resistive RAM (RRAM)** systems for in-memory computing applications.

---

##  Key Features

-  High-performance simulation engine written in C
-  Custom netlist parser and circuit representation
-  Numerical solvers for circuit analysis
-  Supports large circuits (IBM) implemented with sparse matrices

### Supported Analyses:
- DC Analysis  
- AC Analysis  
- Transient Analysis  
- Operating Point (OP) Calculation  

### Supported Components:
- Resistors  
- Capacitors  
- Inductors  
- Diodes  
- Independent Voltage Sources  
- Independent Current Sources  
- **Memristors (native implementation)**  

---

## Memristor Modeling

The simulator implements a native memristor model that eliminates the need for external behavioral definitions.

Users only need to specify model parameters, making simulation significantly more efficient compared to conventional approaches.

---

## Project Structure

```bash
memristor-spice-simulator/
│
├── src/                  # Core simulator source code
│   ├── parser/           # Netlist parser
│   ├── solver/           # Numerical solvers
│   ├── components/       # Circuit elements implementation
│   └── analyses/         # DC, AC, transient, OP analysis
│
├── tests/                # Example netlists
│
├── results/              # Simulation outputs
│
├── Makefile
└── README.md
