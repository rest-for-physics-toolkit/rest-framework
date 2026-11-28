// Version selector of the published documentation (GitHub Pages, branch gh-pages): the site has a folder per
// version (master, vX.Y.Z) and, at its root, this script and versions.json (the list of versions). The docs
// workflow loads it in each page; it replaces the version number of the title area with a list of the
// versions, opening the same page in the version chosen (or its main page if the page does not exist there)
(function() {
const root = new URL('.', document.currentScript.src);

function addSelector(versions) {
    const relative = window.location.href.substring(root.href.length);
    const current = relative.split('/')[0];
    const page = relative.substring(current.length + 1) || 'index.html';
    const place = document.getElementById('projectnumber');
    if (!place || !versions.includes(current)) return;

    const select = document.createElement('select');
    for (const version of versions) {
        const option = document.createElement('option');
        option.value = version;
        option.textContent = version;
        option.selected = version === current;
        select.appendChild(option);
    }
    select.addEventListener('change', function() {
        const target = new URL(select.value + '/' + page, root);
        fetch(target, {method: 'HEAD'})
            .then(function(response) {
                window.location.href = response.ok ? target.href : new URL(select.value + '/', root).href;
            })
            .catch(function() {
                window.location.href = new URL(select.value + '/', root).href;
            });
    });
    place.textContent = ' ';
    place.appendChild(select);
}

const versions = fetch(new URL('versions.json', root)).then(function(response) {
    return response.json();
});
document.addEventListener('DOMContentLoaded', function() {
    versions.then(addSelector).catch(function() {});
});
})();
