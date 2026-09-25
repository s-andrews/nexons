// Embedded in both QC reports so their controls travel with the HTML file.
(() => {
    const ids = ['innerflexchart', 'startflexchart', 'endflexchart'];
    Chart.register({
        id: 'flex-visible-range',
        afterDataLimits(chart, {scale}) {
            if (scale.id !== 'y' || !ids.includes(chart.canvas.id)) return;
            const left = chart.options.scales.x.min ?? -Infinity;
            const right = chart.options.scales.x.max ?? Infinity;
            let peak = 0;
            chart.data.datasets.forEach((dataset, index) => {
                if (!chart.isDatasetVisible(index)) return;
                let previous = null;
                for (const point of dataset.data) {
                    if (!Number.isFinite(point.x) || !Number.isFinite(point.y)) {
                        previous = null;
                        continue;
                    }
                    if (point.x >= left && point.x <= right) peak = Math.max(peak, point.y);
                    // Include lines crossing the viewport edges, even between sparse points.
                    if (previous && point.x > previous.x) {
                        for (const edge of [left, right]) {
                            if (previous.x < edge && edge < point.x) {
                                const y = previous.y + (point.y - previous.y) *
                                    (edge - previous.x) / (point.x - previous.x);
                                peak = Math.max(peak, y);
                            }
                        }
                    }
                    previous = point;
                }
            });
            scale.min = 0;
            scale.max = peak > 0 ? peak * 1.05 : 1;
        }
    });

    for (const id of ids) {
        const chart = Chart.getChart(id);
        const canvas = chart.canvas;
        chart.options.scales.y.min = 0;
        chart.update('none');
        const fullMin = chart.scales.x.min;
        const fullMax = chart.scales.x.max;
        const fullWidth = fullMax - fullMin;
        let drag = null;

        function setRange(min, width) {
            width = Math.min(fullWidth, Math.max(Math.min(1, fullWidth), width));
            min = Math.max(fullMin, Math.min(fullMax - width, min));
            chart.options.scales.x.min = min;
            chart.options.scales.x.max = min + width;
            chart.update('none');
        }
        function position(event) {
            const rect = canvas.getBoundingClientRect();
            return {
                x: (event.clientX - rect.left) * chart.width / rect.width,
                y: (event.clientY - rect.top) * chart.height / rect.height
            };
        }
        function inside({x, y}) {
            const area = chart.chartArea;
            return x >= area.left && x <= area.right && y >= area.top && y <= area.bottom;
        }
        canvas.style.cursor = 'grab';
        canvas.addEventListener('wheel', event => {
            const point = position(event);
            if (!inside(point) || drag) return;
            event.preventDefault();
            const axis = chart.scales.x;
            const fraction = (point.x - axis.left) / axis.width;
            const delta = event.deltaY * (event.deltaMode === 1 ? 16 : event.deltaMode === 2 ? chart.height : 1);
            const width = Math.min(fullWidth, Math.max(Math.min(1, fullWidth),
                (axis.max - axis.min) * Math.exp(Math.max(-1, Math.min(1, delta * 0.002)))));
            setRange(axis.getValueForPixel(point.x) - fraction * width, width);
        }, {passive: false});
        canvas.addEventListener('pointerdown', event => {
            if (event.pointerType !== 'mouse' || event.button !== 0 || !inside(position(event))) return;
            const rect = canvas.getBoundingClientRect();
            drag = {id: event.pointerId, x: event.clientX, min: chart.scales.x.min,
                width: chart.scales.x.max - chart.scales.x.min,
                pixels: chart.scales.x.width * rect.width / chart.width};
            canvas.setPointerCapture(event.pointerId);
            canvas.style.cursor = 'grabbing';
            event.preventDefault();
        });
        canvas.addEventListener('pointermove', event => {
            if (!drag || event.pointerId !== drag.id) return;
            setRange(drag.min - (event.clientX - drag.x) * drag.width / drag.pixels, drag.width);
        });
        function endDrag(event) {
            if (!drag || event.pointerId !== drag.id) return;
            drag = null;
            canvas.style.cursor = 'grab';
            if (canvas.hasPointerCapture(event.pointerId)) canvas.releasePointerCapture(event.pointerId);
        }
        canvas.addEventListener('pointerup', endDrag);
        canvas.addEventListener('pointercancel', endDrag);
        canvas.addEventListener('lostpointercapture', endDrag);

        const controls = document.createElement('div');
        controls.className = 'my-2';
        const reset = document.createElement('button');
        reset.type = 'button';
        reset.className = 'btn btn-sm btn-outline-secondary me-2';
        reset.textContent = 'Reset view';
        reset.setAttribute('aria-controls', id);
        reset.addEventListener('click', () => {
            delete chart.options.scales.x.min;
            delete chart.options.scales.x.max;
            chart.update('none');
        });
        controls.appendChild(reset);
        const hint = document.createElement('small');
        hint.className = 'text-muted';
        hint.textContent = 'Scroll over the plot to zoom; drag left or right to pan.';
        controls.appendChild(hint);
        canvas.parentElement.insertAdjacentElement('afterend', controls);
    }
})();
