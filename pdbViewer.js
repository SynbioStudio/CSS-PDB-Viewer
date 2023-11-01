const EXPECTED_DISTANCE = 25;
    const MAX_DIAMETER = 20; // Maximum size of the residue
    const MIN_DIAMETER = 10; // Minimum size of the residue

    async function fetchPDBData() {
        
const pdbElement = document.querySelector('[data-url-pdb]');
const pdbFileName = pdbElement ? pdbElement.getAttribute('data-url-pdb') : '2B3P.pdb';
const response = await fetch('https://files.rcsb.org/download/' + pdbFileName);

        const text = await response.text();
        return text.split('\n')
            .filter(line => line.startsWith('ATOM') && line.slice(12, 16).trim() === 'CA')
            .map(line => ({
                x: parseFloat(line.slice(30, 38)),
                y: parseFloat(line.slice(38, 46)),
                z: parseFloat(line.slice(46, 54)),
                aminoAcid: line.slice(17, 20).trim(),
                residueNumber: parseInt(line.slice(22, 26).trim())
            }));
    }

    function normalizeCoordinates(residues) {
        const margin = 20;
        const centerX = residues.reduce((acc, r) => acc + r.x, 0) / residues.length;
        const centerY = residues.reduce((acc, r) => acc + r.y, 0) / residues.length;
        const centerZ = residues.reduce((acc, r) => acc + r.z, 0) / residues.length;
        const width = Math.max(...residues.map(r => r.x)) - Math.min(...residues.map(r => r.x));
        const height = Math.max(...residues.map(r => r.y)) - Math.min(...residues.map(r => r.y));
        const depth = Math.max(...residues.map(r => r.z)) - Math.min(...residues.map(r => r.z));

        return residues.map(residue => ({
            x: margin + (residue.x - centerX + width / 2) / width * (400 - 2 * margin),
            y: margin + (residue.y - centerY + height / 2) / height * (400 - 2 * margin),
            z: margin + (residue.z - centerZ + depth / 2) / depth * (400 - 2 * margin),
            aminoAcid: residue.aminoAcid,
            residueNumber: residue.residueNumber
        }));
    }

    function rotate(coordinates, angleX, angleY) {
        let sinX = Math.sin(angleX);
        let cosX = Math.cos(angleX);
        let sinY = Math.sin(angleY);
        let cosY = Math.cos(angleY);

        // Calculate center of mass
        let centerX = coordinates.reduce((acc, r) => acc + r.x, 0) / coordinates.length;
        let centerY = coordinates.reduce((acc, r) => acc + r.y, 0) / coordinates.length;
        let centerZ = coordinates.reduce((acc, r) => acc + r.z, 0) / coordinates.length;

        return coordinates.map(coord => {
            // Translate to origin
            let dx = coord.x - centerX;
            let dy = coord.y - centerY;
            let dz = coord.z - centerZ;

            // Perform rotation
            let ry1 = cosX * dy - sinX * dz;
            let rz1 = sinX * dy + cosX * dz;
            let rx2 = cosY * dx + sinY * rz1;
            let rz2 = -sinY * dx + cosY * rz1;

            // Translate back to the original position
            return { 
                x: rx2 + centerX, 
                y: ry1 + centerY, 
                z: rz2 + centerZ, 
                aminoAcid: coord.aminoAcid, 
                residueNumber: coord.residueNumber 
            };
        });
    }

    function renderPDB(residues) {
        const viewer = document.getElementById('pdbViewer');
        while (viewer.firstChild) {
            viewer.removeChild(viewer.firstChild);
        }

        residues.forEach((residue, index) => {
            const perspectiveFactor = 1 - (residue.z / 400);
            const diameter = MIN_DIAMETER + (MAX_DIAMETER - MIN_DIAMETER) * perspectiveFactor;
            const residueElement = document.createElement('div');
            residueElement.classList.add('residue');
            residueElement.style.width = `${diameter}px`;
            residueElement.style.height = `${diameter}px`;
            residueElement.style.left = `${residue.x}px`;
            residueElement.style.top = `${residue.y}px`;
            residueElement.setAttribute('data-info', `${residue.aminoAcid} (${residue.residueNumber})`);
            residueElement.style.zIndex = `${Math.round(1000 - residue.z)}`;
            viewer.appendChild(residueElement);

            if (index < residues.length - 1) {
                const nextResidue = residues[index + 1];
                const deltaX = nextResidue.x - residue.x;
                const deltaY = nextResidue.y - residue.y;
                const angle = Math.atan2(deltaY, deltaX) * (180 / Math.PI);
                const distance = Math.sqrt(deltaX ** 2 + deltaY ** 2) - (diameter / 2) - (diameter / 2);
                const bondElement = document.createElement('div');
                bondElement.classList.add('bond');
                bondElement.style.width = `${distance}px`;
                bondElement.style.left = `${residue.x + (diameter / 2) * Math.cos(angle * Math.PI / 180)}px`;
                bondElement.style.top = `${residue.y + (diameter / 2) * Math.sin(angle * Math.PI / 180)}px`;
                bondElement.style.transform = `rotate(${angle}deg)`;
                bondElement.style.zIndex = `${Math.round((residue.z + nextResidue.z) / 2 * 1000)}`;
                if (distance > EXPECTED_DISTANCE) {
                    bondElement.classList.add('dotted');
                }
                viewer.appendChild(bondElement);
            }

    // Add mouseenter event to dynamically adjust z-index
    residueElement.addEventListener('mouseenter', function() {
        this.style.zIndex = '999999999';  // Set to a very high value

    });
        });
    }

    let isDragging = false;
    let lastMouseX = null;
    let lastMouseY = null;
    const rotationSpeed = 0.05;
    let currentRotationX = 0;
    let currentRotationY = 0;

    document.getElementById('pdbViewer').addEventListener('mousedown', (e) => {
        isDragging = true;
        lastMouseX = e.clientX;
        lastMouseY = e.clientY;
    });

    document.addEventListener('mouseup', () => {
        isDragging = false;
    });

    document.getElementById('pdbViewer').addEventListener('mousemove', (e) => {
        if (!isDragging) return;
        const dx = e.clientX - lastMouseX;
        const dy = e.clientY - lastMouseY;
        currentRotationX += dy * rotationSpeed;
        currentRotationY += dx * rotationSpeed;
        const rotatedResidues = rotate(normalizedResidues, currentRotationX, currentRotationY);
        renderPDB(rotatedResidues);
        lastMouseX = e.clientX;
        lastMouseY = e.clientY;
    });

    let normalizedResidues;
    fetchPDBData().then(residues => {
        normalizedResidues = normalizeCoordinates(residues);
        const initialRotationX = Math.PI / 4;
        const initialRotationY = Math.PI / 4;
        const rotatedResidues = rotate(normalizedResidues, initialRotationX, initialRotationY);
        renderPDB(rotatedResidues);
    });