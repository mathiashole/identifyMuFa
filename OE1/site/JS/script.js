  let stage = new NGL.Stage("viewport", { backgroundColor: "white" });
  let currentComponent;
  let currentPerspective = "perspective";
  let autoSpin = false;  // Estado inicial: auto rotación apagada

//   function cargarModelo(file) {
//     stage.removeAllComponents();
//     stage.loadFile(file).then(function (component) {
//       currentComponent = component;
//       actualizarRepresentacion();
//       stage.autoView();
//       // Al cargar un modelo, respetamos el estado actual de autoSpin
//       stage.setSpin(autoSpin);
//       actualizarBotonAutoSpin();
//     });
//   }
function cargarModelo(file) {
  stage.removeAllComponents();
  stage.loadFile(file).then(function(component) {
    currentComponent = component;
    actualizarRepresentacion();
    stage.autoView();
    stage.setSpin(autoSpin);
    actualizarBotonAutoSpin();
    
    // Obtener y mostrar la secuencia
    mostrarSecuencia(component);
  });
}

function mostrarSecuencia(component) {
  const sequenceContainer = document.getElementById('sequence-viewer');
  sequenceContainer.innerHTML = '';
  
  // Obtener la estructura
  const structure = component.structure;
  
  // Crear contenedor principal
  const container = document.createElement('div');
  container.className = 'sequence-container';
  
  // Obtener las cadenas (chains)
  const chains = structure.getChainnames();
  
  chains.forEach(chainId => {
    // Crear sección para cada cadena
    const chainDiv = document.createElement('div');
    chainDiv.className = 'sequence-chain';
    
    // Añadir etiqueta de la cadena
    const label = document.createElement('div');
    label.className = 'sequence-label';
    label.textContent = `Cadena ${chainId}`;
    chainDiv.appendChild(label);
    
    // Obtener residuos para esta cadena
    const residues = [];
    structure.eachResidue(function(residue) {
      if (residue.chainname === chainId) {
        residues.push(residue);
      }
    });
    
    // Mostrar secuencia
    const seqDiv = document.createElement('div');
    residues.forEach(residue => {
      const resElem = document.createElement('span');
      resElem.className = 'residue';
      resElem.textContent = residue.resname;
      resElem.title = `Residuo ${residue.resno}: ${residue.resname}`;
      
      // Resaltar al hacer hover
      resElem.addEventListener('mouseover', function() {
        // Resaltar residuo en la visualización 3D
        currentComponent.addRepresentation('ball+stick', {
          sele: `:${residue.resno} and .${chainId}`,
          color: 'red'
        });
        
        // Resaltar en la secuencia
        resElem.classList.add('highlight');
      });
      
      resElem.addEventListener('mouseout', function() {
        // Quitar resaltado
        currentComponent.removeRepresentation(currentComponent.reprList.length - 1);
        resElem.classList.remove('highlight');
      });
      
      seqDiv.appendChild(resElem);
    });
    
    chainDiv.appendChild(seqDiv);
    container.appendChild(chainDiv);
  });
  
  sequenceContainer.appendChild(container);
}


  function actualizarRepresentacion() {
    if (!currentComponent) return;
    currentComponent.removeAllRepresentations();
    const rep = document.getElementById("rep-select").value;
    const color = document.getElementById("color-select").value;
    currentComponent.addRepresentation(rep, { colorScheme: color });
    currentComponent.autoView();
  }

  function toggleAutoSpin() {
    autoSpin = !autoSpin;
    stage.setSpin(autoSpin);
    actualizarBotonAutoSpin();
  }

  function actualizarBotonAutoSpin() {
    const btn = document.getElementById("autoSpinBtn");
    btn.textContent = autoSpin ? "🔄 Auto-rotación: ON" : "🔄 Auto-rotación: OFF";
  }

  document.getElementById("model-select").addEventListener("change", function () {
    cargarModelo(this.value);
  });

  document.getElementById("color-select").addEventListener("change", actualizarRepresentacion);
  document.getElementById("rep-select").addEventListener("change", actualizarRepresentacion);

  function descargarImagen() {
    stage.makeImage({ factor: 2 }).then(function (blob) {
      NGL.download(blob, "estructura.png");
    });
  }

  function pantallaCompleta() {
    const elem = document.getElementById("viewport");
    if (elem.requestFullscreen) {
      elem.requestFullscreen();
    } else if (elem.webkitRequestFullscreen) {
      elem.webkitRequestFullscreen();
    } else if (elem.msRequestFullscreen) {
      elem.msRequestFullscreen();
    }
  }

  function cambiarPerspectiva() {
    currentPerspective = (currentPerspective === "perspective") ? "orthographic" : "perspective";
    stage.setParameters({ cameraType: currentPerspective });
  }

  // Cargar modelo inicial
  cargarModelo(document.getElementById("model-select").value);