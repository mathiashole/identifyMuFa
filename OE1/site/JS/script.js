  let stage = new NGL.Stage("viewport", { backgroundColor: "white" });
  let currentComponent;
  let currentPerspective = "perspective";
  let autoSpin = false;  // Estado inicial: auto rotación apagada

  function cargarModelo(file) {
    stage.removeAllComponents();
    stage.loadFile(file).then(function (component) {
      currentComponent = component;
      actualizarRepresentacion();
      stage.autoView();
      // Al cargar un modelo, respetamos el estado actual de autoSpin
      stage.setSpin(autoSpin);
      actualizarBotonAutoSpin();
    });
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